! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE PRINT_EXTRA(LOOKUP,D1,NS_PTS,K)

      USE lookup_addresses

      IMPLICIT NONE

! Purpose: To print out extra data from D1
!
! Subroutine Arguments:

      INTEGER                                                           &
     & LOOKUP(64,*)                                                     &
                          ! IN  :lookup table
     &,NS_PTS                                                           &
                          ! IN  :number of extra data values printed
     &,K                  ! IN  :field number

      REAL                                                              &
     & D1(*)              ! IN  :D1 array

      INTEGER                                                           &
     & TOT_VALUES                                                       &
     &,IEXTRAW                                                          &
     &,ADDR                                                             &
     &,CODE                                                             &
     &,DATA_VALUES                                                      &
     &,MAX_PRINT                                                        &
     &,EW_PTS                                                           &
                          ! No of points across to print
     &,EW_PRINT                                                         &
     &,START_OFF                                                        &
     &,N_VECT                                                           &
                          ! No of vector in extra data
     &,MAX_INTHD                                                        &
     &,MAX_DATA                                                         &
                          ! Max possible size of extra data
     &,I,J,L                                                            &
                          !Loop counter
     &,ICODE

      PARAMETER( EW_PTS=5,                                              &
     &           MAX_INTHD=30)


      INTEGER IND   ! Local index


      INTEGER                                                           &
     & X_INTHD(MAX_INTHD)

      REAL                                                              &
     & DATA(LOOKUP(LBEXT,K))                                            &
     &,DATA_LINE(EW_PTS)

      CHARACTER(LEN=12)                                                      &
     & DASH(EW_PTS)       !Stores dashed lines

      CHARACTER(LEN=1)                                                       &
     & BLANK

      PARAMETER(BLANK=' ')

      CHARACTER(LEN=12)                                                      &
     & CMESSAGE


! Extract information from extra data
! Put the vector header of each extra data vector
! into array X_INTHD
! Put the first NS_PTS values of each extra data vector
! into array X_DATA

! A extra data vector comprises of:
!
!  header value1 value2 value3 value4 ....
!
! |______|______|______|______|______|_____|
!
!
! The extra data vector(s) are attached to the end of
! each set of field data.  For a detail description,
! see documentation F3
!

      IEXTRAW=LOOKUP(LBEXT,K)
      TOT_VALUES=LOOKUP(LBLREC,K)
      ADDR=TOT_VALUES-IEXTRAW+1
      L=0
      IND=0
      START_OFF=1
      ! While our address is short of the final record
      DO WHILE (ADDR <  TOT_VALUES)
        L=L+1
        ! Decode the extra data vector descriptor
        CODE=TRANSFER(D1(ADDR), CODE)
        X_INTHD(L)=CODE
! DEPENDS ON: check_extra
        CALL CHECK_EXTRA(CODE,DATA_VALUES,ICODE,CMESSAGE)
        ICODE=0
        ! The max number of points to print out - don't exceed
        ! the number of elements in each vector
        MAX_PRINT=MIN(DATA_VALUES,NS_PTS)

        ! Fill in actual values
        I = 0
        DO J=ADDR+1,ADDR+DATA_VALUES
          I = I + 1
          IF (I <= MAX_PRINT) THEN
             IND=IND+1
             DATA(IND)=D1(J)
          ENDIF
        ENDDO
        ADDR=ADDR+DATA_VALUES+1   ! INCREMENT ADDRESS
                                  ! by no of data elements + 1
                                  ! for the vector descriptor.
        START_OFF=START_OFF+MAX_PRINT
      ENDDO
      N_VECT=L

! Set up print format

      DO J=1,EW_PTS
        DASH(J)='------------'
      ENDDO

      EW_PRINT=MIN(N_VECT,EW_PTS)

! Print out extra data

! Print 1st row - integer extra data header

      WRITE(6,'(14X,A25,I3,A8)')                                        &
     &   'Extra Data Values (first ',MAX_PRINT,' values)'
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':'',9(I10,2X))')                     &
     &     K,(X_INTHD(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

! Print remaining rows - extra data values
! For now, just printing out max of 5 vectors

      DO START_OFF=1,MAX_PRINT
        DO I=1,EW_PRINT
          DATA_LINE(I)=DATA(START_OFF+(I-1)*MAX_PRINT)
        ENDDO
        WRITE(6,'(14X,A1,9(F9.3,3X))')                                  &
     &  BLANK,(DATA_LINE(J),J=1,EW_PRINT)
      ENDDO

      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' '')')

      IF(N_VECT >  EW_PTS) THEN
        WRITE(6,*)'NUMBER OF EXTRA DATA VECTORS: ',                     &

     &  N_VECT
        WRITE(6,*)'ONLY THE FIRST 5 WERE PRINTED OUT.'
        WRITE(6,*)' '
      ENDIF

      RETURN
      END SUBROUTINE PRINT_EXTRA
