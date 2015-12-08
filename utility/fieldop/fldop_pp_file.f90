! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!Subroutine interface:
      SUBROUTINE FLDOP_PP_FILE(PPFIELD,LENBUF,NUM_WORDS,                &
     &RMDI,COMP_ACCRCY,                                                 &
     &PPHORIZ_OUT,UNITPP,DATA_ADD,N_COLS_OUT,N_ROWS_OUT,PACKING,        &
     &PACKING_TYPE,LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,ENTRY_NO,             &
     &ICODE,CMESSAGE)
      USE IO

      USE io_configuration_mod, ONLY  : &
          io_field_padding

      IMPLICIT NONE
!
! Description: To output a field to a PP_FILE
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     &  ICODE                                                           &
                           !   RETURN CODE FROM ROUTINE
     &, LENBUF                                                          &
                           !   LENGTH OFF PP BUFFER
     &, UNITPP                                                          &
                           !   OUTPUT PP UNIT NUMBER
     &, LEN_IO             !NOT USED, BUT NEEDED FOR BUFFOUT CALL

      INTEGER                                                           &
     &  N_ROWS_OUT                                                      &
                      !  PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT
     &, N_COLS_OUT                                                      &
                      !   PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT
     &, NUM_OUT                                                         &
                      !   NUMBER OF COMPRESSED (32 BIT) WORDS
     &, DATA_ADD                                                        &
                      !
     &, ENTRY_NO                                                        &
                      !
     &, COMP_ACCRCY                                                     &
                      !   PACKING ACCURACY IN POWER OF 2
     &, PPHORIZ_OUT                                                     &
                      !   SIZE OF OUTPUT FIELD
     &, NUM_WORDS                                                       &
                      !   NUMBER OF 64 BIT WORDS WORTH OF DATA
     &, PACKING_TYPE                                                    &
     &, LEN1_LOOKUP                                                     &
     &, LEN2_LOOKUP

      REAL                                                              &
     & RMDI                   !IN     Missing data indicator

      LOGICAL                                                           &
     &  PACKING            !IN OVERALL Packing switch (T if pckng reqd)

!   Array  arguments with intent(in):
      INTEGER                                                           &
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)

      REAL                                                              &
     & PPFIELD(PPHORIZ_OUT)   !INOUT ARRAY TO STORE PPDATA

!   Scalar arguments with intent(out):
      CHARACTER(LEN=80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE

!   Array arguments with intent(out):
      REAL                                                              &
     & BUFOUT(LENBUF)         !OUTPUT PP BUFFER (ROUNDED UP)

! Local scalars:
      INTEGER                                                           &
     & LENGTH_FULLWRD                                                   &
                     !     LENGTH IN BITS OF FULLWORD VAR
     &,LEN_BUF_WORDS                                                    &
                     !     NUM_WORDS ROUNDED BY 512
     &,POS

      INTEGER                                                           &
     &  JJ            !     Local counter

      REAL                                                              &
     &  IX            !     RETURN VALUE FROM UNIT COMMAND

!- End of header


      LENGTH_FULLWRD=64   !   LENGTH IN BITS OF FULLWORD VAR
      IF(PACKING_TYPE == 1 .AND. PACKING)THEN

! DEPENDS ON: coex
        CALL COEX(PPFIELD,PPHORIZ_OUT,BUFOUT,LENBUF,N_COLS_OUT,         &
     &  N_ROWS_OUT,NUM_OUT,COMP_ACCRCY,.TRUE.,RMDI,LENGTH_FULLWRD,      &
     &  ICODE,CMESSAGE)

        NUM_WORDS=(NUM_OUT+1)/2 ! Round up to the nearest 64 Bit CRAY Wd
        LEN_BUF_WORDS=((NUM_WORDS+IO_FIELD_PADDING-1)/IO_FIELD_PADDING)*    &
     &                                                IO_FIELD_PADDING

      ELSE IF(PACKING_TYPE == 4 .AND. PACKING) THEN
! DEPENDS ON: runlen_encode
          CALL RUNLEN_ENCODE(PPFIELD,PPHORIZ_OUT,BUFOUT,PPHORIZ_OUT,    &
     &                     NUM_OUT,RMDI,ICODE,CMESSAGE)
          NUM_WORDS=NUM_OUT
          LEN_BUF_WORDS=((NUM_WORDS+IO_FIELD_PADDING-1)/IO_FIELD_PADDING)*  &
     &    IO_FIELD_PADDING
      ELSE  ! No packing required.

        DO JJ=1,PPHORIZ_OUT
          BUFOUT(JJ) = PPFIELD(JJ)
        END DO

        NUM_WORDS=PPHORIZ_OUT
        LEN_BUF_WORDS=((NUM_WORDS+IO_FIELD_PADDING-1)/IO_FIELD_PADDING)*    &
     &                                 IO_FIELD_PADDING

      ENDIF

      ! Update lookup header data lengths and addressing for
      ! wgdos packed data in fieldsfile.
      LOOKUP(15,ENTRY_NO) = NUM_WORDS
      LOOKUP(30,ENTRY_NO) = LEN_BUF_WORDS
      IF (ENTRY_NO  ==  1) THEN
        LOOKUP(29,ENTRY_NO) = DATA_ADD
      ELSE
        LOOKUP(29,ENTRY_NO) = LOOKUP(29,ENTRY_NO-1)                     &
     &                        + LOOKUP(30,ENTRY_NO-1)
      ENDIF
      LOOKUP(40,ENTRY_NO) = LOOKUP(29,ENTRY_NO)
      ! Set position in output file to buffer out lookup header info.
      POS = LOOKUP(40,ENTRY_NO)

      DO JJ=NUM_WORDS+1,LEN_BUF_WORDS
        BUFOUT(JJ)= 0.0
      ENDDO

      CALL SETPOS(UNITPP,POS,ICODE)
      CALL BUFFOUT(UNITPP,BUFOUT(1:),LEN_BUF_WORDS,LEN_IO,IX)

      RETURN
      END SUBROUTINE FLDOP_PP_FILE
