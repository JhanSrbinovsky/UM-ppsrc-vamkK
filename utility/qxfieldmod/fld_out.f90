! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!LL  SUBROUTINE FLD_OUT------------------------------------------
!LL
!LL  REPLACES THE OUTPUT FROM STASH EITHER ON TO A PP FILE OR
!LL  BACK TO THE MAIN ARRAY D1
!LL
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL  VERSION 1, DATED 12/09/89
!LL
!LL  SYSTEM TASK: CONTROL PART OF C4
!LL
!LL  PURPOSE:   TO PROCESS DIAGNOSTICS CONTROLLED BY STASH
!LL
!LL
!LLEND-------------------------------------------------------------

!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs


!
!*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE FLD_OUT                                                &
     &          (ICODE,CMESSAGE,BUFOUT,LENBUF,LEN_BUF_WORDS,NUM_WORDS,  &
     &           UNITPP,LEN1_LOOKUP,PP_LEN2_LOOKUP,IPPLOOK,RPPLOOK,     &
     &           ILABEL,RLABEL,IWL,DATA_ADDR)
      USE IO
      IMPLICIT NONE

      CHARACTER(LEN=*) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
!

      INTEGER                                                           &
     &  ICODE                                                           &
                           !IN    RETURN CODE FROM ROUTINE
     &, LEN1_LOOKUP                                                     &
                           !IN    FIRST DIMENSION OF LOOKUP TABLE
     &, PP_LEN2_LOOKUP                                                  &
                           !IN    SECND DIMENSION OF LOOKUP TABLE
     &, LENBUF                                                          &
                           !IN     LENGTH OFF PP BUFFER
     &, UNITPP                                                          &
                           !IN     OUTPUT PP UNIT NUMBER
     &, LEN_BUF_WORDS                                                   &
                           !IN
     &, NUM_WORDS          !IN
!
      INTEGER                                                           &
     &  JJ            !IN    ITEM NUMBER
      INTEGER                                                           &
     &  IPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP)                             &
                                            !IN INTEGER LOOKUP TABLE
     &, ILABEL(45)                                                      &
                      ! INTEGER PART OF LOOKUP
     &, IWL                                                             &
                      !IN    Address of the PP LOOKUP Table
     &, DATA_ADDR     !IN    Address of start of data
!
      REAL                                                              &
     &  BUFOUT(LENBUF)                                                  &
                               !OUTPUT PP BUFFER (ROUNDED UP)
     &, RPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP)                             &
                                            !IN REAL LOOKUP TABLE
     &, RLABEL(19)    ! REAL PART OF LOOKUP

!*---------------------------------------------------------------------

!*L  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
!*---------------------------------------------------------------------
!     EQUIVALENCE(IPPLOOK,RPPLOOK)
!
!*------------------------------------------------------------------
!L  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
!L---------------------------------------------------------------------
!----------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES
      INTEGER                                                           &
     &  ADDR                                                            &
                      !
     &, IWA                                                             &
                      !     RECORD NUMBER
     &, IX                                                              &
                      !     RETURN VALUE FROM UNIT COMMAND
     &, LEN_IO                                                          &
                      !
     &, II                                                              &
                      !     COUNTER
     &, I                                                               &
                      !     COUNTER
     &,IERR           ! Error return from SETPOS

      real                                                              &
     &  A_IO          !

      INTEGER                                                           &
     &  LRESID                                                          &
                      !
     &, ICURRLL                                                         &
                      !
     &, IPAST                                                           &
                      !
     &, IPROJ         !     M08 PROJECTION NUMBER

      LOGICAL                                                           &
     &  FIRST              !
      DATA FIRST/.TRUE./

      LOGICAL :: found
!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
      FIRST=.TRUE.
      ICURRLL=0

  103 FORMAT(//,32X,' ARRAY FROM START OF PPOUT  ',//,32(10F8.0/))
      LRESID=LEN_BUF_WORDS-NUM_WORDS
      DO JJ=NUM_WORDS+1,LRESID
      BUFOUT(JJ)= 0.0
      END DO
!
      found=.TRUE.
      IF(FIRST) THEN
        found=.FALSE.
        DO  JJ=1,PP_LEN2_LOOKUP
          IF(IPPLOOK(1,JJ) <  0) THEN  ! Search for last entry
            ICURRLL=JJ
            IF(JJ == 1) THEN
              IWA=((IWL+511)/512)*512+PP_LEN2_LOOKUP*LEN1_LOOKUP
              write(6,*) 'Start data',iwa,data_addr,iwa-1
              IWA=DATA_ADDR
              IWA=IWA-1
            ELSE
              IWA= IPPLOOK(29,JJ-1)+IPPLOOK(30,JJ-1) !ADDR+LGTH
            ENDIF
            found=.true.
            EXIT
          ENDIF
        END DO
        IF (.NOT.found) THEN
          ICODE=1
          CMESSAGE="FROM PPOUT CANNOT FIND SUITABLE ENTRY IN LOOKUP"
          ! FIXME : Should there be an ereport here?
        END IF
      ELSE
          IPAST=ICURRLL-1
        WRITE(7,105) IPAST
 105    FORMAT('  FROM PPOUT AND FIRST IS FALSE IPAST=',I8)
          IWA=IPPLOOK(29,IPAST) + IPPLOOK(30,IPAST) ! ADDR + LENGTH
        WRITE(7,106) IWA
 106    FORMAT('  FROM PPOUT AND FIRST IS FALSE IWA=',I8)
      ENDIF

      IF (found) THEN
!
!     update lookup for this field
        DO I=1,45
          IPPLOOK(I,ICURRLL) = ILABEL(I)
        ENDDO
        DO I=1,19
          RPPLOOK(I+45,ICURRLL) = RLABEL(I)
        ENDDO
        IPPLOOK(29,ICURRLL)=IWA
        IPPLOOK(30,ICURRLL)=LEN_BUF_WORDS
        IPPLOOK(40,ICURRLL)=IWA
        !
        CALL SETPOS(UNITPP,IWA,IERR)
        CALL BUFFOUT(unitpp,bufout,LEN_buf_words,LEN_IO,A_IO)
100     FORMAT(//,32X,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10F8.0/))
101     FORMAT(//,32X,'   IPPLOOK AT END OF  PPOUT   ',//,32(16I5/))
102     FORMAT('     IWA  LEN_BUF_WORDS ',2I12)
      END IF
      RETURN
      END SUBROUTINE FLD_OUT
