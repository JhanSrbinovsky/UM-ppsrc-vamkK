! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   SUBROUTINE COEX,COEX2,CMPS,XPND,INSTIN,EXTRIN -----------------
!
!   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!
!
!   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!   STANDARD B, VERSION 2, DATED 18/01/90
!
!  Logical component number: S72
!
!   SYSTEM TASK: P7
!  -------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: STASH
      SUBROUTINE XPND(IX,ICOMP,FIELD,APREC,IBIT,BASE,NOP,RMDI)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

!     Subroutine arguments
      INTEGER IX,NOP,ICOMP(NOP+1),IBIT
      REAL FIELD(IX),APREC,BASE,RMDI
      REAL ATEMP(IX)

!     Local variables
      LOGICAL OBTMIN,OBTMIS,OBTZER,OBTMAP
      INTEGER :: I,JWORD,JBIT,JJ,NPOINT,II,IEND
      INTEGER :: IMAP(IX),IMIS(IX),IZERO(IX),IMIN(IX)
      INTEGER :: full_word, iword, idiv, jword_end, jbit_end
      INTEGER IGATH(IX)




!     FUNCTIONS
      INTEGER ishift

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!
!         IX   =  LENGTH OF FIELD
!      ICOMP   =  DATA FOR EXPANSION
!      FIELD   =  FIELD OF EXPANDED DATA
!      APREC   =  PRECISION
!       IBIT   =  NUMBER OF BITS INTO WHICH DATA PACKED
!       BASE   =  REFERENCE (MIN) VALUE FOR PACKED DATA
!        NOP   =  SIZE OF COMP
!       RMDI   =  MISSING DATA INDICATOR
!
!
!     INITIALISE VARIABLES/TEMP ARRAYS
!
      IF (lhook) CALL dr_hook('XPND',zhook_in,zhook_handle)
      OBTMAP   = .FALSE.
      OBTMIS   = .FALSE.
      OBTMIN   = .FALSE.
      OBTZER   = .FALSE.
      DO I = 1,IX
          IMAP(I)  = 1
          IMIS(I)  = 0
          IMIN(I)  = 0
          IZERO(I) = 0
      END DO
!
!     CHECK IF BITMAP USED FOR ZERO VALUES
!
      IF (IBIT >= 128) THEN
          OBTZER = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT-128
      ELSE
          OBTMIN = .FALSE.
      ENDIF
!
!     CHECK IF BITMAP USED FOR MINIMUM VALUES (NOT CURRENTLY USED)
!
      IF (IBIT >= 64) THEN
          OBTMIN = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT - 64
      ELSE
          OBTMIN = .FALSE.
      ENDIF
!
!     CHECK IF BITMAP USED FOR MISSING DATA
!
      IF (IBIT >= 32) THEN
          OBTMIS = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT - 32
      ELSE
          OBTMIS = .FALSE.
      ENDIF
!
!     SET START POSITION IN ICOMP
!
      JWORD = 1
      JBIT  = 63
!
!     EXTRACT BITMAPS
!
      IF (OBTMIS) THEN
!
!         EXTRACT MISSING DATA BITMAP
!
          DO I=1,IX
              IF (BTEST(ICOMP(JWORD),JBIT)) THEN
                  IMIS(I) = 1
                  IMAP(I) = 0
              ELSE
                  IMIS(I) = 0
              ENDIF
              IF(JBIT >  0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
          END DO
      ENDIF
      IF (OBTMIN) THEN
!
!         EXTRACT MINIMUM VALUE BITMAP (NOT USED AT PRESENT)
!
          DO I=1,IX
              IF(BTEST(ICOMP(JWORD),JBIT)) THEN
                  IMIN(I) = 1
                  IMAP(I) = 0
              ELSE
                  IMIN(I) = 0
              ENDIF
              IF(JBIT >  0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
          END DO
      ENDIF
      IF (OBTZER) THEN
!
!         EXTRACT ZERO VALUE BITMAP
!
!         For faster execution we first check if all 64 bits in a word
!         are zero before we test the bits separately.
!
          IEND = MIN(IX,JBIT+1)
          DO I=1,IEND
              IF(BTEST(ICOMP(JWORD),JBIT-I+1)) THEN
                  IZERO(I)= 0
              ELSE
                  IZERO(I)= 1
                  IMAP(I) = 0
              ENDIF
          END DO
          JBIT = JBIT - IEND
          IF(JBIT <  0) THEN
              JBIT    = 63
              JWORD   = JWORD + 1
          ENDIF
          DO II=IEND+1,IX,64
              IF (ICOMP(JWORD) == 0) THEN
                  IZERO(II:MIN(IX,II+64-1)) = 1
                  IMAP (II:MIN(IX,II+64-1)) = 0
              ELSE
                  DO I=II,MIN(IX,II+64-1)
                      IF(BTEST(ICOMP(JWORD),JBIT-I+II)) THEN
                          IZERO(I)= 0
                      ELSE
                          IZERO(I)= 1
                          IMAP(I) = 0
                      ENDIF
                  END DO
              ENDIF
              JBIT = JBIT - ( MIN(IX,II+64-1) - II + 1 )
              IF(JBIT <  0) THEN
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
          END DO
      ENDIF
!
!     IF BIT MAP USED FIND NUMBER OF POINTS STORED
!     AND RESET POINTERS TO BEGINNING OF 32 BIT BOUNDARY
!
      IF (OBTMAP) THEN
          NPOINT = 0
          DO I=1,IX
              IF (IMAP(I) == 1) NPOINT = NPOINT + 1
          END DO
          IF (JBIT /= 63) THEN
              IF (JBIT >= 31) THEN
                  JBIT  = 31
              ELSE
                  JBIT  = 63
                  JWORD = JWORD + 1
              ENDIF
          ENDIF
      ELSE
          NPOINT = IX
      ENDIF
      IF (IBIT >  0) THEN
!
!         UNPACK SCALED VALUES TO TEMP ARRAY
!
          jword_end = jword + (NPOINT*IBIT)/64
          jbit_end  = jbit  - MOD(NPOINT*IBIT,64)
          IF (MOD(IBIT,64) == 0) THEN
              IDIV = 64
          ELSE IF (MOD(IBIT,32) == 0) THEN
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
!CDIR NODEP
          iword = jword
          DO I = II+1,NPOINT,64/IDIV
!
! Get the significant bits from the first word into the
! top of 'full_word'
!
! DEPENDS ON: ishift
              full_word=ishift(icomp(iword), 63-jbit)
!
! Shift the second word to generate spaces for the part from the
! first word, and then combine the two parts in 'full_word'
!
              full_word=ior(full_word,                                  &
! DEPENDS ON: ishift
     &         ishift(icomp(iword+1), -(jbit+1)))
!
! Now shift 'full_word' down so that the required integer is in
! the bottom most bits
!
! DEPENDS ON: ishift
              full_word=ishift(full_word, -(64-ibit))
!
!         ADD DIFFERENCES TO MINIMUM VALUE AND UNSCALE
!
              ATEMP(I)=full_word*APREC+BASE
              iword = iword + IBIT/IDIV
          ENDDO
!
          JBIT = JBIT - IBIT
          IF (JBIT <  0) THEN
              JWORD = JWORD + 1
              JBIT  = JBIT + 64
          ENDIF

          ENDDO
!
! Calculate the position in the source array
!
          jword = jword_end
          jbit  = jbit_end
          IF (JBIT <  0) THEN
              JWORD = JWORD + 1
              JBIT  = JBIT + 64
          END IF
!
!
!         MOVE INTO UNPACKED ARRAY FIELD
!
! FIRST GET GATHER INDEX
!
          JJ  =0
          DO I=1,IX
            IF (IMAP(I) == 1) THEN
              JJ       =JJ+1
              IGATH(JJ)=I
            END IF
          END DO
!
          DO I=1,IX
          FIELD(I)=0.
          END DO

          DO I=1,JJ
          FIELD(IGATH(I)) = ATEMP(I)
          END DO
!
!         IF MINIMUMS BIT MAPPED FILL ZEROS IN FIELD WITH BASE
!
          IF (OBTMIN) THEN
              DO I=1,IX
                  IF(IMIN(I) == 1) FIELD(I) = BASE
              END DO
          ENDIF
!
!         IF MISSING DATA BIT MAPPED FILL ZEROS IN FIELD WITH MDI
!
          IF (OBTMIS) THEN
              DO I=1,IX
                  IF(IMIS(I) == 1) FIELD(I) = RMDI
              END DO
          ENDIF

      ELSEIF (IBIT == 0) THEN

!
!         ALL POINTS IN ROW HAVE SAME VALUE E.G. POLE ROW
!
          DO I=1,IX
              FIELD(I)=BASE
          END DO
!
!         IF MISSING DATA BIT MAPPED FILL ZEROS IN FIELD WITH MDI
!
          IF (OBTMIS) THEN
              DO I=1,IX
                  IF(IMIS(I) == 1) FIELD(I) = RMDI
              END DO
          ENDIF
      ENDIF
!
      IF (lhook) CALL dr_hook('XPND',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE XPND
