! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: This routine acts as an interface between the model and
!  GRIB format output routines.
!
!  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!  Logical components covered: ...
!
!  Project task: ...
!
!  -------------------------------------------------------------------
!  Interface and arguments: ------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
      SUBROUTINE GRIB_FILE(LEN1_LOOKUP,LEN2_LOOKUP,LOOKUP,RLOOKUP,IENT, &
     &                     FIELD,PPHORIZ_OUT,LENBUF,NUM_CRAY_WORDS,     &
     &                     UNITPP,IWA,GRIB_PACKING,ICODE,CMESSAGE)
      USE IO
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE lookup_addresses

      IMPLICIT NONE

      INTEGER                                                           &
     &     LEN1_LOOKUP                                                  &
                       !  IN   first dimension of LOOKUP
     &    ,LEN2_LOOKUP                                                  &
                       !  IN   second dimension of LOOKUP
     &    ,LENBUF                                                       &
                       !  IN   No of points in output field
     &    ,IENT                                                         &
                       !  IN   level indicator for processing LOOKUP.
     &    ,IWA                                                          &
                       !  IN   Record number
     &    ,PPHORIZ_OUT                                                  &
                       !  IN
     &    ,UNITPP                                                       &
                       !  IN   Output PP unit number
     &    ,GRIB_PACKING                                                 &
                        !  IN  Packing profile for grib
     &    ,LEN_FIELD                                                    &
     &    ,ICODE                                                        &
                          !  OUT  Return code
     &    ,NUM_CRAY_WORDS                                               &
                          !  OUT  Number of cray words output in grib
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Integer lookup headers

      REAL                                                              &
     &     FIELD(PPHORIZ_OUT)                                           &
                               ! IN   Unpacked output array
     &    ,RLOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! REAL lookup headers
      CHARACTER                                                         &
     &     CMESSAGE*(*)     ! OUT  Will contain any error messages
!
! LOCAL VARIABLES
!
      INTEGER                                                           &
     &     ILABEL(45)                                                   &
                            ! Integer part of LOOKUP for level IENT
     &    ,LEN_IO                                                       &
     &    ,IX
      REAL  :: Rerror
      REAL                                                              &
     &     RLABEL(19)                                                   &
                            ! Real part of LOOKUP for level IENT
     &    ,WORK_ARRAY(LENBUF)                                           &
                              ! GRIB packed output array
     &    ,BUFOUT(LENBUF)   ! Output PP BUFFER

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER :: i,j ! Loop indices

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!L
!L 1. Fill arrays ILABEL and RLABEL
!L
      IF (lhook) CALL dr_hook('GRIB_FILE',zhook_in,zhook_handle)
      DO J=1,45
        ILABEL(J)=LOOKUP(J,IENT)
      ENDDO
      DO J=1,19
        RLABEL(J)=RLOOKUP(J+45,IENT)
      ENDDO
!L
!L 2. Convert data to GRIB code
!L
! DEPENDS ON: pp2grib
      CALL PP2GRIB(FIELD,WORK_ARRAY,LENBUF,NUM_CRAY_WORDS,GRIB_PACKING, &
     &             ILABEL,RLABEL,ICODE,CMESSAGE)
      IF (ICODE /= 0) THEN
        IF (lhook) CALL dr_hook('GRIB_FILE',zhook_out,zhook_handle)
        RETURN
      ENDIF
!     WRITE(6,*) NUM_CRAY_WORDS,LENBUF
!     write(6,*) (ilabel(j),j=1,45)
!     write(6,*) (rlabel(j),j=1,19)
!L
!L 3. Put coded data into BUFOUT for output
!L
      DO I=1,NUM_CRAY_WORDS
        BUFOUT(I)=WORK_ARRAY(I)
      ENDDO
      DO I=NUM_CRAY_WORDS+1,LENBUF
        BUFOUT(I)=0.0
      ENDDO
!L
!L 4. Update lookup for this field
!L
      DO J=1,45
        LOOKUP(J,IENT)=ILABEL(J)
      ENDDO
      DO J=1,19
        RLOOKUP(J+45,IENT)=RLABEL(J)
      ENDDO
      LOOKUP(LBLREC,IENT)=NUM_CRAY_WORDS
      LOOKUP(LBEGIN,IENT)=IWA
      LOOKUP(LBNREC,IENT)=NUM_CRAY_WORDS
      LOOKUP(DATA_TYPE,IENT)=1
      LOOKUP(NADDR,IENT)=IWA
!L
!L 5. Output BUFOUT
!L
      CALL SETPOS(UNITPP,IWA,ICODE)
      CALL BUFFOUT(UNITPP,BUFOUT,NUM_CRAY_WORDS,LEN_IO,Rerror)
      IF (lhook) CALL dr_hook('GRIB_FILE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GRIB_FILE
