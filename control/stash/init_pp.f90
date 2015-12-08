! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! 
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: STASH
!
!  Purpose: Initialises direct access PP files at the start of
!           the run.  NB: Sequential PP files need no initialisation.
!
!  Logical components covered: D401
!
!  External documentation: On-line UM document C61 - Zonal mean
!                          calculations.



      SUBROUTINE INIT_PP ( FTN_UNIT,FILE_TYPE_LETTER,                   &
     &                     LEN1_LOOKUP,PP_LEN2_LOOKUP,FIXHD,            &
     &                     INTHD,REALHD,LEVDEPC,ROWDEPC,COLDEPC,        &
     &                     LEN_INTHD,                                   &
     &                     LEN_REALHD,LEN1_LEVDEPC,LEN2_LEVDEPC,        &
     &                     LEN1_ROWDEPC,LEN2_ROWDEPC,                   &
     &                     LEN1_COLDEPC,LEN2_COLDEPC,                   &
     &                     PP_LEN_INTHD, PP_LEN_REALHD,                 &
     &                     ICODE,CMESSAGE)


      USE Model_file                                                 

      USE IOS, ONLY : IOS_Set_Header
      USE io_configuration_mod, ONLY :   &
           io_data_alignment
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE IO, ONLY                   :   &
          ioDiskSynchronise,             &
          buffout
      USE UM_ParVars

      IMPLICIT NONE
!
      CHARACTER(LEN=1)                                                       &
     &    FILE_TYPE_LETTER    ! IN  - File type (p-PP, b-bndry)
      INTEGER                                                           &
     &    FTN_UNIT                                                      &
                              ! IN  - Fortran unit number
     &,   LEN1_LOOKUP                                                   &
                              ! IN  - Size of PP header
     &,   PP_LEN2_LOOKUP                                                &
                              ! IN  - Max allowable fields
     &,   LEN_INTHD                                                     &
                              ! IN    LENGTH OF INTEGER CONSTANTS
     &,   LEN_REALHD                                                    &
                              ! IN    LENGTH OF REAL CONSTANTS
     &,   LEN1_LEVDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of lev depndt
     &,   LEN2_LEVDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of lev depndt
     &,   LEN1_ROWDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of row depndt
     &,   LEN2_ROWDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of row depndt 
     &,   LEN1_COLDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of col depndt
     &,   LEN2_COLDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of col depndt
     &,   ICODE                                                         &
                              ! OUT - Error exit code
     &,   PP_LEN_INTHD                                                  &
                              ! IN - Length of PP FILE integer header
     &,   PP_LEN_REALHD       ! IN - Length of PP FILE real header
!
!
      INTEGER                                                           &
     &    FIXHD(*)                                                      &
                                    ! IN    ARRAY OF FIXED CONSTANTS
     &,   INTHD(LEN_INTHD)
                                    ! IN    ARRAY OF integer CONSTANTS
      REAL                                                              &
     &    LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)                            &
                                              ! IN LEV DEP CONSTANTS
     &,   ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC)                            &
                                              ! IN ROW DEP CONSTANTS
     &,   COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC)
                                              ! IN COL DEP CONSTANT
      INTEGER                                                           &
     &    PP_INTHD(PP_LEN_INTHD)
                                    ! OUT   ARRAY of integer constants
      REAL                                                              &
     &    PP_LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)                         &
                                                 ! OUT Level dep cts
     &,   PP_ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC)                         &
                                                 ! OUT Row dep cts
     &,   PP_COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC)  ! OUT Col dep cts


      INTEGER, POINTER :: PP_FIXHD(:)





!
      REAL                                                              &
     &    REALHD(LEN_REALHD)                                            &
                                    ! IN    ARRAY OF REAL CONSTANTS
     &,   PP_REALHD(PP_LEN_REALHD)  ! OUT   ARRAY OF REAL CONSTANTS
!
      CHARACTER(LEN=80)                                                      &
     &    CMESSAGE            ! OUT - Error message
!
!*----------------------------------------------------------------------
!
!
!  Local variables
!

      INTEGER, POINTER   :: IPPLOOK(:,:)
      INTEGER, PARAMETER :: current_io_pe=0

      INTEGER            :: DUMMY
      INTEGER            :: STEP



!
!dir$ cache_align pp_fixhd, pp_inthd, pp_realhd, pp_levdepc, 
!     pp_rowdepc, pp_coldepc, ipplook
      INTEGER                                                           &
     &       II,JJ,IWA,IX,LEN_IO,START_BLOCK  !
      REAL A_IO

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
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

      IF (lhook) CALL dr_hook('INIT_PP',zhook_in,zhook_handle)
!L----------------------------------------------------------------------
!L 1. Reserve space
!L
       NULLIFY(IPPLOOK)
       STEP = 1
       CALL initLookups(IPPLOOK, FTN_UNIT, LEN1_LOOKUP,                 &
     &                   PP_LEN2_LOOKUP, DUMMY, STEP)

!L----------------------------------------------------------------------
!L 1.1 Set up FIXED header record for the PP FILE
!L
! Attach fixed length header
      CALL initHeader(ftn_unit,FixedHeader)
      PP_FIXHD=>attachHeader(ftn_unit,FixedHeader)
      DO II=1,SIZE(pp_fixhd)
        PP_FIXHD(II)=FIXHD(II)
      END DO
      IF (FILE_TYPE_LETTER == 'p' .OR.                                  &
     &    FILE_TYPE_LETTER == 'c') THEN
        PP_FIXHD(5)=3
      ELSEIF (FILE_TYPE_LETTER == 'b') THEN
        PP_FIXHD(5)=5
      ELSE
        ICODE=100
        CMESSAGE='INIT_PP  : Unknown output file type letter'
        IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
        RETURN
      ENDIF
      PP_FIXHD(101)=PP_LEN_INTHD
      PP_FIXHD(105)=PP_FIXHD(100)+PP_FIXHD(101)
      PP_FIXHD(106)=PP_LEN_REALHD
      PP_FIXHD(110)=PP_FIXHD(105)+PP_FIXHD(106)
      PP_FIXHD(111)=LEN1_LEVDEPC
      PP_FIXHD(112)=LEN2_LEVDEPC
      PP_FIXHD(115)=0
      PP_FIXHD(116)=imdi
      PP_FIXHD(117)=imdi
      PP_FIXHD(120)=0
      PP_FIXHD(121)=imdi
      PP_FIXHD(122)=imdi
      PP_FIXHD(125)=0
      PP_FIXHD(126)=imdi
      PP_FIXHD(127)=imdi
      PP_FIXHD(130)=0
      PP_FIXHD(131)=imdi
      PP_FIXHD(135)=0
      PP_FIXHD(136)=imdi
      PP_FIXHD(140)=0
      PP_FIXHD(141)=imdi
      PP_FIXHD(142)=0
      PP_FIXHD(143)=imdi
      PP_FIXHD(144)=0
      PP_FIXHD(145)=imdi
      PP_FIXHD(150)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
      IF (LEN2_ROWDEPC > 0) THEN
        PP_FIXHD(115)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
        PP_FIXHD(116)=LEN1_ROWDEPC 
        PP_FIXHD(117)=LEN2_ROWDEPC
        PP_FIXHD(150)=PP_FIXHD(115)+ PP_FIXHD(116)*PP_FIXHD(117)  
      END IF 
      IF (LEN2_COLDEPC > 0) THEN 
        PP_FIXHD(120)=PP_FIXHD(115)+ PP_FIXHD(116)*PP_FIXHD(117)
        PP_FIXHD(121)=LEN1_COLDEPC
        PP_FIXHD(122)=LEN2_COLDEPC
        PP_FIXHD(150)=PP_FIXHD(120)+ PP_FIXHD(121)*PP_FIXHD(122) 
      END IF
      PP_FIXHD(151)=LEN1_LOOKUP
      PP_FIXHD(152)=PP_LEN2_LOOKUP
      pp_fixhd(160)=                                                    &
                         ! make sure the data starts on a sector bndry
     & ((pp_fixhd(150)+pp_len2_lookup*len1_lookup-1+io_data_alignment-1)/ &
     & io_data_alignment)*io_data_alignment+1

      CALL SetHeader(ftn_unit,FixedHeader)

!L----------------------------------------------------------------------
!L 1.2 Set up INTEGER constants record for the PP FILE
!L
      IF(PP_FIXHD(5) <= 2) THEN !  set all values initially to MDI
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=INTHD(21)
        ENDDO
      ELSE
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=IMDI
        ENDDO
      ENDIF

      PP_INTHD(6)=INTHD(6)
      PP_INTHD(7)=INTHD(7)
      PP_INTHD(8)=INTHD(8)
      PP_INTHD(9)=INTHD(9)
      PP_INTHD(10)=INTHD(10)
      PP_INTHD(12)=INTHD(12)
      PP_INTHD(13)=INTHD(13)
      PP_INTHD(17)=INTHD(17)
      PP_INTHD(24)=INTHD(24)
      PP_INTHD(25)=INTHD(25)
      PP_INTHD(28)=INTHD(28)
!L----------------------------------------------------------------------
!L 1.3 Set up REAL constants record for the PP FILE
!L
      DO II = 1, PP_LEN_REALHD
        PP_REALHD(II) = RMDI   ! Set all values to RMDI initially
      END DO
      PP_REALHD(1)=REALHD(1)
      PP_REALHD(2)=REALHD(2)
      PP_REALHD(3)=REALHD(3)
      PP_REALHD(4)=REALHD(4)
! Set to RMDI for VR      
      IF (LEN2_ROWDEPC > 0 .AND. LEN2_COLDEPC > 0) THEN
        PP_REALHD(1) = RMDI
        PP_REALHD(2) = RMDI     
        PP_REALHD(3) = RMDI    
        PP_REALHD(4) = RMDI       
      ENDIF      
      PP_REALHD(5)=REALHD(5)
      PP_REALHD(6)=REALHD(6)

      PP_REALHD(16)=REALHD(16)
      PP_REALHD(17)=REALHD(17)




!L----------------------------------------------------------------------
!L 1.4 Set up LEVEL/ROW/COL DEPENDANT constants record for the PP FILE
!L
      DO II=1,LEN1_LEVDEPC*LEN2_LEVDEPC
        PP_LEVDEPC(II)=LEVDEPC(II)
      END DO
      
      DO II=1,LEN1_ROWDEPC*LEN2_ROWDEPC
         PP_ROWDEPC(II)=ROWDEPC(II)
      END DO
      
      DO II=1,LEN1_COLDEPC*LEN2_COLDEPC
         PP_COLDEPC(II)=COLDEPC(II)
      END DO
       
!L----------------------------------------------------------------------
!L 2.1 BUFFER OUT Header Records starting with the FIXED LENGTH
!L

      CALL BUFFOUT(FTN_UNIT,PP_FIXHD,SIZE(pp_fixhd),LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= SIZE(pp_fixhd)) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('bufferout of fixed length header',A_IO,LEN_IO, &
     &                    SIZE(pp_fixhd))
           CMESSAGE='INIT_PP:I/O error'
           ICODE=1
           IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
           RETURN
        ENDIF
      START_BLOCK=SIZE(pp_fixhd)+1
!L----------------------------------------------------------------------
!L 2.2 BUFFER OUT Integer Constants
!L

      IF(FIXHD(100) >  0) THEN  ! Any integer constants to output ?

! Check for error in file pointers

!        WRITE(6,*)  'START_BLOCK FIXHD(100)'
!        WRITE(6,*)   START_BLOCK
!        WRITE(6,*)   FIXHD(100)
!        WRITE(6,*)   FTN_UNIT
         IF(FIXHD(100) /= START_BLOCK) THEN  ! Check start address
! DEPENDS ON: poserror
            CALL POSERROR('integer constants',START_BLOCK,100,          &
     &      PP_FIXHD(100))
            CMESSAGE='INIT_PP:  Addressing conflict'
            ICODE=2
            IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
            RETURN
         END IF


         CALL BUFFOUT (FTN_UNIT,PP_INTHD(1:),PP_FIXHD(101),LEN_IO,A_IO)

! Check for I/O errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_FIXHD(101)) THEN
! DEPENDS ON: ioerror
            CALL IOERROR('buffer out of integer constants',A_IO,LEN_IO  &
     &     ,PP_FIXHD(101))
            CMESSAGE='INIT_PP: I/O Error'
            ICODE=3
            IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
            RETURN
         END IF

         START_BLOCK=START_BLOCK+PP_FIXHD(101)

      END IF

!L----------------------------------------------------------------------
!L 2.3 BUFFER OUT Real Constants
!L

      IF(PP_FIXHD(105) >  0) THEN   ! Any real constants to output ?

! Check for error in file pointers

        IF(PP_FIXHD(105) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,PP_FIXHD(105))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=4
          IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
          RETURN
        END IF


        CALL BUFFOUT(FTN_UNIT,PP_REALHD(1:),PP_FIXHD(106),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= PP_FIXHD(106)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of real constants',A_IO,LEN_IO       &
     &                 ,PP_FIXHD(106))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=5
          IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
          RETURN
        END IF

        START_BLOCK=START_BLOCK+PP_FIXHD(106)

      END IF

!L----------------------------------------------------------------------
!L 2.4.1 BUFFER OUT Level Dependant Constants.
!L

      IF(PP_FIXHD(112) >  0) THEN ! Any level dependant constants ?

! Check for error in file pointers

         IF(PP_FIXHD(110) /= START_BLOCK) THEN
! DEPENDS ON: poserror
            CALL POSERROR('real constants',START_BLOCK,100,             &
     &                     PP_FIXHD(110))
            CMESSAGE='INIT_PP: Addressing conflict'
            ICODE=6
            IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
            RETURN
         END IF


         CALL BUFFOUT (FTN_UNIT,PP_LEVDEPC(1:)                           &
     &              ,PP_FIXHD(111)*PP_FIXHD(112),LEN_IO,A_IO)

! Check for I/O errors

         IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(111)*PP_FIXHD(112)      &
     &        ))THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer out of lev dep constants',A_IO,LEN_IO   &
     &            ,PP_FIXHD(111))
           CMESSAGE='INIT_PP: I/O Error'
           ICODE=7
           IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
           RETURN
         END IF

         START_BLOCK=START_BLOCK+ PP_FIXHD(111)*PP_FIXHD(112)

      END IF
!L----------------------------------------------------------------------
!L 2.4.2 BUFFER OUT Row Dependant Constants.
!L

      IF(PP_FIXHD(115) >  0) THEN ! Any row dependant constants ?

! Check for error in file pointers

        IF(PP_FIXHD(115) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,               &
     &                   PP_FIXHD(115))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=6
          IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
          RETURN
        END IF


        CALL BUFFOUT (FTN_UNIT,PP_ROWDEPC(1:),                           &
     &                PP_FIXHD(116)*PP_FIXHD(117),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(116)*PP_FIXHD(117)       &
     &       ))THEN

! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of row dep constants',A_IO,LEN_IO,   &
     &                  PP_FIXHD(116))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=7
          IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
          RETURN
        END IF

        START_BLOCK=START_BLOCK+ PP_FIXHD(116)*PP_FIXHD(117)

      END IF
!L----------------------------------------------------------------------
!L 2.4.3 BUFFER OUT Col Dependant Constants.
!L

      IF(PP_FIXHD(120) >  0) THEN ! Any col dependant constants ?

! Check for error in file pointers

        IF(PP_FIXHD(120) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,               &
     &                   PP_FIXHD(120))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=6
          IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
          RETURN
        END IF


        CALL BUFFOUT (FTN_UNIT,PP_COLDEPC(1:),                           &
     &                PP_FIXHD(121)*PP_FIXHD(122),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(121)*PP_FIXHD(122)       &
     &       ))THEN

! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of col dep constants',A_IO,LEN_IO ,  &
     &                  PP_FIXHD(121))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=7
          IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
          RETURN
        END IF

        START_BLOCK=START_BLOCK+ PP_FIXHD(121)*PP_FIXHD(122)

      END IF      
!L----------------------------------------------------------------------
!L 2.5 BUFFER OUT Lookup Table
!L
!     IWA= 0
!     CALL SETPOS(FTN_UNIT,3,IWA,ICODE)
           IF(PP_FIXHD(152) >  0) THEN

! Check for error in file pointers

             IF(PP_FIXHD(150) /= START_BLOCK) THEN
! DEPENDS ON: poserror
               CALL POSERROR('lookup table',START_BLOCK,100,            &
     &              PP_FIXHD(150))
               CMESSAGE='INIT_PP: Addressing conflict'
               ICODE=8
               IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
               RETURN
             END IF


      CALL BUFFOUT (FTN_UNIT,                                           &
     &              IPPLOOK,LEN1_LOOKUP*PP_LEN2_LOOKUP,LEN_IO,A_IO)

!
! Check for I/O errors

            IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(151)*PP_FIXHD(152))) &
     &          THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of PP LOOKUP TABLE ',A_IO,LEN_IO &
     &            ,PP_FIXHD(152))
              CMESSAGE='INIT_PP: I/O Error'
              ICODE=9
              IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
              RETURN
            END IF

            ! Force last write to then file to commit to disk
            ! to avoid problems with continuation runs following 
            ! hard failures.
            CALL ioDiskSynchronise(ftn_unit)

            START_BLOCK=START_BLOCK+(PP_FIXHD(151)*PP_FIXHD(152))

!
! If we are the current I/O PE, we need to update our copy
! of the LOOKUP Table disk address
!
      STEP = 2
      CALL initLookups(IPPLOOK, FTN_UNIT, DUMMY, DUMMY,                &
     &                  PP_FIXHD(150) - 1, STEP)
      NULLIFY(PP_FIXHD)


          END IF
      IF (lhook) CALL dr_hook('INIT_PP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INIT_PP
