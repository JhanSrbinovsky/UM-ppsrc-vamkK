! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE COEX,COEX2,CMPS,XPND,INSTIN,EXTRIN -----------------
!LL
!LL   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!LL
!LL
!LL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL   STANDARD B, VERSION 2, DATED 18/01/90
!LL
!LL  Logical component number: S72
!LL
!LL   SYSTEM TASK: P7
!LL   OWNED BY P J SMITH
!LL
!LLEND-------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: STASH
      SUBROUTINE COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,          &
     &                 ICODE,CMESSAGE)

      USE ereport_mod, ONLY : ereport
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

!     Subroutine arguments
      INTEGER N,ICOMP(N),M,IX,IY,NUM,ISC
      INTEGER ICODE
      CHARACTER CMESSAGE*80
      REAL FIELD(M),RMDI
      LOGICAL OCO

!     Local variables
      INTEGER JJ,IST,ICX,JCX,NOB,IER2
      INTEGER IC(IY),ICB(IY),NOP(IY),IBIT(IY),ISTART(IY)
      INTEGER JCOMP(IX,IY),IERR1(IY),IERR2(IY),IERR3(IY)
      Integer :: Ip
      Integer :: Ierr
      Real :: Acc
      Real :: Aprec
      Real :: Base(Iy)

      Character (Len=*), Parameter :: routinename='coex2'
!     External functions used
      Integer :: IEEE2IBM
      Integer :: IBM2IEEE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!     CRAY VERSION C.1.1  16/11/90  P J SMITH
!     CRAY ROUTINE     64 BIT WORDS - CRAY FULLWORD
!
!                     OCO=.TRUE.                 OCO=.FALSE.
!
!      FIELD   =  FIELD FOR COMPRESSING     RETURNED EXPANDED DATA
!          M   =  SIZE OF FIELD             SIZE OF FIELD
!      ICOMP   =  RETURNED COMPRESSED DATA  FIELD FOR EXPANSION
!          N   =  SIZE OF COMP                 -------
!         IX   =  X DIMENSION OF FIELD      X DIMENSION OF FIELD
!         IY   =  Y DIMENSION OF FIELD      Y DIMENSION OF FIELD
!        NUM   =  TOTAL NO. OF COMPRESSED      -------
!                 (32 BIT) WORDS RETURNED
!        ISC   =  ACCURACY IN POWER OF 2    ACCURACY IN POWER OF 2
!        OCO   =  .TRUE. FOR COMPRESSION    .FALSE. FOR EXPANSION
!        IRC   =  RETURN CODE FOR EACH ROW  RETURN CODE FOR EACH ROW
!
!     INITIALISE TEMP. VARIABLES/ARRAYS
!
      IF (lhook) CALL dr_hook('COEX2',zhook_in,zhook_handle)
      IF (.NOT. OCO) THEN
        DO JJ=1,IY
          NOP(JJ)=0
          IBIT(JJ)=0
        END DO
        DO JJ=1,IY
          BASE(JJ)=0.0
          IBIT(JJ)=0.0
        END DO
      END IF

      DO JJ=1,IY
          DO JCX=1,IX
             JCOMP(JCX,JJ) = 0
          END DO
      END DO
      IF (OCO) THEN
!
!         COMPRESSION OF DATA
!
          IC(1) = 1
          ACC   = ISC
          APREC = 2.**ACC
!
!         PUT PACKED DATA FOR EACH ROW INTO TEMP ARRAY JCOMP
!

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                 &
!$OMP& PRIVATE(jj, ip, cmessage)                                  &
!$OMP& SHARED(ix, iy, field, jcomp, nop, aprec, ibit, base, rmdi, &
!$OMP&        ierr1, ierr2, ierr3)                                &
!$OMP& FIRSTPRIVATE(icode)
          DO JJ = 1,IY
!             IP      = POSITION IN INPUT ARRAY
              IP      = (JJ-1)*IX + 1

! DEPENDS ON: cmps
              CALL CMPS(IX,FIELD(IP),JCOMP(2,JJ),NOP(JJ),APREC,         &
     &         IBIT(JJ),BASE(JJ),RMDI,ICODE,CMESSAGE)

            IF (ICODE  /=  0) THEN
              CALL ereport(routinename, icode, cmessage)
            END IF
              

!
!         ADD BASE VALUE, NUMBER OF BITS, AND LENGTH
!
              IERR1(JJ)=IEEE2IBM(3,1, JCOMP(1,JJ),0, BASE(JJ),1,64,32)
              IERR2(JJ)=IEEE2IBM(2,1,JCOMP(1,JJ),32,IBIT(JJ),1,64,16)
              IERR3(JJ)=IEEE2IBM(2,1,JCOMP(1,JJ),48,NOP(JJ),1,64,16)
          END DO
!$OMP END PARALLEL DO
!
!         CHECK ROW HEADER AND SET RETURN CODES
!
!         CALCULATE POSITIONS IN OUTPUT ARRAY FOR PACKED ROWS
!         (FIRST PACKED ROW STARTS AT WORD 1; BIT 31)
!
          IC(1)     = 2
          ICB(1)    = -1
          ISTART(1) = 5
          DO JJ = 2,IY
              IF (MOD(NOP(JJ-1),2) == 1) THEN
                  IC(JJ ) = IC(JJ-1) + NOP(JJ-1)/2 + 1
                  ICB(JJ) = -ICB(JJ-1)
                  IF (ICB(JJ) >  0) IC(JJ) = IC(JJ) + 1
              ELSE
                  IC(JJ)  = IC(JJ-1) + (NOP(JJ-1)+1)/2 + 1
                  ICB(JJ) = ICB(JJ-1)
              ENDIF
              ISTART(JJ)  = 5
              IF(ICB(JJ) == 1) ISTART(JJ) = 1
          END DO
!
!         MOVE TEMP. ARRAY INTO OUTPUT ARRAY
!
          DO JJ=1,IY
              NOB  = NOP(JJ)*4 + 8
! CHECK IF PACKED FIELD GREATER THAN UN PACKED FIELD
!             IF(NOB >  IX*8)IRC(JJ)=IRC(JJ)+32
              IST  = ISTART(JJ)
              ICX  = IC(JJ)
              CALL MOVEBYTES(JCOMP(1,JJ),1,NOB,ICOMP(ICX),IST)
          END DO
!
!         INSERT TOTAL LENGTH OF THIS FIELD
          NUM = IC(IY)*2 + NOP(IY)
          IF (ICB(IY) <  0) NUM = NUM + 1
          IER2=IEEE2IBM(2,1,ICOMP(1),0,NUM,1,64,32)
!
!         END OF COMPRESSION SECTION
!
      ELSE
!
!         EXPANSION SECTION
!
          ACC   = ISC
          APREC = 2.**ACC
          ICX   = 2
          JCX   = -1

          DO JJ = 1,IY
!
!             MOVE PACKED ROWS INTO TEMP ARRAYS
!
              IF (JCX <  0) THEN
!
!                 EXTRACT BASE, NO. BITS, NO 32 BIT WORDS
!
                  IERR=IBM2IEEE(3,1,ICOMP(ICX),32,BASE(JJ),1,64,32)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX+1),0,IBIT(JJ),1,64,16)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX+1),16,NOP(JJ),1,64,16)
!                 SAVE START POSITION OF ROW
                  IC(JJ)     = ICX
                  ISTART(JJ) = 5
              ELSE
                  IERR=IBM2IEEE(3,1,ICOMP(ICX),0,BASE(JJ),1,64,32)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX),32,IBIT(JJ),1,64,16)
                  IERR=IBM2IEEE(2,1,ICOMP(ICX),48,NOP(JJ),1,64,16)
!                 SAVE START POSITION OF ROW
                  IC(JJ)     = ICX
                  ISTART(JJ) = 1
              END IF
!
!             CALCULATE START OF NEXT ROW
!
              IF (MOD(NOP(JJ),2) == 1) THEN
                  ICX   = ICX + NOP(JJ)/2 + 1
                  JCX   = -JCX
                  IF (JCX >  0) ICX = ICX + 1
              ELSE
                  ICX   = ICX + (NOP(JJ)+1)/2 + 1
              END IF
          END DO
!
!         MOVE EACH PACKED ROW INTO TEMP ARRAY JCOMP
!
          DO JJ = 1,IY
              ICX  = IC(JJ)
              IST  = ISTART(JJ)
              NOB  = NOP(JJ)*4 + 8
              CALL MOVEBYTES(ICOMP(ICX),IST,NOB,JCOMP(1,JJ),1)
          END DO
!
!         CALCULATE START OF EACH ROW IN FIELD
!
          ICX = 1
          DO JJ = 1,IY
              IC(JJ) = ICX
              ICX    = ICX + IX
          END DO
!
!         UNPACK DATA INTO FIELD
!
          DO JJ=1,IY
              ICX    = IC(JJ)
! DEPENDS ON: xpnd
              CALL XPND(IX,JCOMP(2,JJ),FIELD(ICX),APREC,IBIT(JJ),       &
     &                                            BASE(JJ),NOP(JJ),RMDI)
          END DO
      END IF
 9990  CONTINUE
      IF (lhook) CALL dr_hook('COEX2',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE COEX2
