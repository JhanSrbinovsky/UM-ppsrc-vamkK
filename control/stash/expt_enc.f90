! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE EXPT_ENC--------------------------------------------
!
!     Given a valid five character experiment RUN_ID, an INTEGER*4
!   code number is generated. The valid experiment name characters
!   are A-Z (uppercase), 0-9. Each letter of the experiment name is
!   stored as a 6 bit number. The last letter is converted to the
!   6 least significant bits of the integer code, the next letter
!   the next 6 lsb's, etc. Hence 30 bits are used in all.
!     The lookup table is capable of holding 64 elements. This
!   number cannot be exceeded if the code number is to remain
!   INTEGER*4. Similarly, the experiment RUN_ID length (5 chars)
!   cannot be exceeded.
!     Subroutine called from PP_HEAD.
!
!   Programming standard:
!
!   Logical components covered:
!
!   Project TASK:
!
!   External documentation:
!
!  -------------------------------------------------------------

!  INTERFACE and ARGUMENTS:--------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE EXPT_ENC(EXPTSTRG                                      &
     &  ,EXPTCODE                                                       &
     &  ,ICODE                                                          &
     &  ,CMESSAGE)
!*-----------------------------------------------------------------


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER       ICODE       !OUT  Return code: successful=0
      CHARACTER(LEN=80)  CMESSAGE    !OUT  Error message if ICODE > 0

      CHARACTER(LEN=5)   EXPTSTRG    !IN   Experiment name string. Length
                                !     must equal parameter STRSIZE
      INTEGER       EXPTCODE    !OUT  Experiment code integer

!     Define local variables

      LOGICAL       TEST
      CHARACTER(LEN=1)   LETTER
      INTEGER       I,J,NEWNUM,LETNUM,NSTRINGS,NBITS,STRSIZE

      PARAMETER(NSTRINGS=36,                                            &
     &  NBITS=6,                                                        &
     &  STRSIZE=5)

      CHARACTER(LEN=1)   USTRINGS    ! Upper case strings
      CHARACTER(LEN=1)   LSTRINGS    ! Lower case strings

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      DIMENSION     USTRINGS(0:NSTRINGS-1)
      DIMENSION     LSTRINGS(0:NSTRINGS-1)

      DATA USTRINGS/'A','B','C','D','E','F','G','H','I','J',            &
     &  'K','L','M','N','O','P','Q','R','S','T',                        &
     &  'U','V','W','X','Y','Z','0','1','2','3',                        &
     &  '4','5','6','7','8','9'/

      DATA LSTRINGS/'a','b','c','d','e','f','g','h','i','j',            &
     &  'k','l','m','n','o','p','q','r','s','t',                        &
     &  'u','v','w','x','y','z','0','1','2','3',                        &
     &  '4','5','6','7','8','9'/

!     Begin main

      IF (lhook) CALL dr_hook('EXPT_ENC',zhook_in,zhook_handle)
      EXPTCODE=0
      LETNUM=STRSIZE

!     Loop over letters in EXPTSTRG
      DO I=0,STRSIZE-1
        TEST=.FALSE.
        READ(EXPTSTRG(LETNUM:LETNUM),"(A1)")LETTER

!       Loop over letters in lookup table USTRINGS/LSTRINGS
        DO J=0,NSTRINGS-1
          IF ((LETTER == USTRINGS(J)).OR.(LETTER == LSTRINGS(J))) THEN
            TEST=.TRUE.
            NEWNUM=J*(2**(I*NBITS))
            EXPTCODE=EXPTCODE+NEWNUM
!           Exit loop as we have found the code for this letter
            EXIT
          ENDIF
        END DO

!       Check experiment name is valid
        IF (.NOT.TEST) THEN
          ICODE=99
          CMESSAGE='EXPT_ENC: Invalid letter in expt name (RUN_ID)'
          GOTO 9999
        ENDIF
        LETNUM=LETNUM-1
      END DO

 9999 CONTINUE
      IF (lhook) CALL dr_hook('EXPT_ENC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE EXPT_ENC
