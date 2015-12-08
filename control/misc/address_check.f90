! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine ADDRESS_CHECK --------------------------------------
!
! Purpose : Check that start addresses of fields read in agree
!           with start addresses set up by UI. Called in INITDUMP
!           if prognostic fields read in from atmos or ocean dumps.
!
! Coding Standard : UM documentation paper no. 3
!
! Documentation : None
!
!----------------------------------------------------------------
!
!*L Arguments
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc

      SUBROUTINE ADDRESS_CHECK (LOOKUP,MPP_DUMP_ADDR,MPP_DUMP_LEN,      &
     &                          LEN1_LOOKUP,LEN2_LOOKUP,                &
     &                          SI,NITEMS,NSECTS,LEN_DATA,              &
     &                          ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Submodel_Mod

      IMPLICIT NONE


      INTEGER                                                           &
     &    LEN1_LOOKUP                                                   &
                              !  1st dimension of lookup table
     &   ,LEN2_LOOKUP                                                   &
                              !  2nd dimension of lookup table
     &   ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                               &
                                           !  Lookup table
     &   ,MPP_DUMP_ADDR(LEN2_LOOKUP)                                    &
                                      ! Addresses of fields as
!                                     ! calculated in READDUMP.
     &   ,MPP_DUMP_LEN(LEN2_LOOKUP)                                     &
                                      ! Lengths of fields as
!                                     ! calculated in READDUMP.

     &   ,LEN_DATA                                                      &
                              !  Expected length of data
     &   ,NITEMS                                                        &
                              !  No of stash items
     &   ,NSECTS                                                        &
                              !  No of stash sections
!  Stash item addresses
     &   ,SI(NITEMS,0:NSECTS,N_INTERNAL_MODEL)                          &
     &   ,ICODE               !  Return code

      CHARACTER(LEN=80)                                                    &
     &          CMESSAGE      !  Error message

!L Dynamic allocation of arrays for F_TYPE
      INTEGER                                                           &
     &    PP_NUM   (LEN2_LOOKUP)                                        &
     &   ,PP_LEN   (LEN2_LOOKUP)                                        &
     &   ,PP_POS   (LEN2_LOOKUP)                                        &
     &   ,PP_LS    (LEN2_LOOKUP)                                        &
     &   ,PP_STASH (LEN2_LOOKUP)                                        &
     &   ,PP_TYPE  (LEN2_LOOKUP)

!L Local array
      INTEGER FIXHD5  !  Dummy variable until removed from F_TYPE

!L Local variables
      INTEGER                                                           &
     &    ADDRESS_STASH                                                 &
     &   ,ADDRESS_LOOKUP                                                &
     &   ,ITEM_CODE                                                     &
     &   ,J                                                             &
     &   ,LEN                                                           &
     &   ,N_TYPES                                                       &
     &   ,SECT_NO                                                       &
     &   ,im_ident                                                      &
                   ! Internal model identifier
     &   ,im_index                                                      &
                   ! Position of int mod id in INTERNAL_MODEL_LIST
     &,OLD_STASH   ! VALUE OF STASH NUMBER ON PREVIOUS ITERATION OF LOOP

      CHARACTER(LEN=80) TITLE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!*---------------------------------------------------------------------
!
!   SET INITIAL VALUE OF PREVIOUS STASH NUMBER
!
      IF (lhook) CALL dr_hook('ADDRESS_CHECK',zhook_in,zhook_handle)
      OLD_STASH = -1
!
!L   Internal Structure

      TITLE = 'Prognostic fields'
! Set fixhd5 to be a dump to allow any checks in f_type to be safe.
      FIXHD5 = 1 
! DEPENDS ON: f_type
      CALL F_TYPE (LOOKUP,LEN2_LOOKUP,PP_NUM,N_TYPES,PP_LEN,            &
     &             PP_STASH,PP_TYPE,PP_POS,PP_LS,FIXHD5,                &
     &TITLE,.FALSE.)

      DO J=1,N_TYPES

!       Get Stash Section no and Item Code
        ITEM_CODE = MOD ( PP_STASH(J),1000)
        SECT_NO   = (PP_STASH(J)-ITEM_CODE)/1000

!       Get im_ident/index for this field
        im_ident = LOOKUP(45,pp_pos(j))
        im_index = INTERNAL_MODEL_INDEX(im_ident)

!       Get lookup and stash start address
        ADDRESS_LOOKUP = MPP_DUMP_ADDR(PP_POS(J))
        ADDRESS_STASH  = SI(ITEM_CODE,SECT_NO,im_index)

!       Check that they match
!
!            CHECK THAT START ADDRESSES AGREE FOR FIRST OCCURRENCE
!            OF A NEW STASH CODE:  FOR FIXED LENGTH FIELDS THERE
!            IS ONLY ONE ENTRY IN THE PP_STASH ARRAY FOR EACH
!            STASH CODE, BUT FOR PACKED FIELDS (EG OCEAN)
!            EACH LEVEL MIGHT HAVE A DIFFERENT LENGTH AND
!            GENERATE A NEW PP_STASH VALUE.
!
         IF (ADDRESS_STASH  /=  ADDRESS_LOOKUP .AND.                    &
     &       OLD_STASH  /=  PP_STASH(J) ) THEN
          CMESSAGE = 'ADDR_CHK : Mis-match in start addresses'
          WRITE (6,*) ' Stash Sect No ',SECT_NO,' Item No ',ITEM_CODE
          WRITE (6,*) ' Start Address in SI           ',ADDRESS_STASH
          WRITE (6,*) ' Start Address in LOOKUP Table ',ADDRESS_LOOKUP
          WRITE (6,*) ' You probably need to RECONFIGURE the start dump'
          ICODE = 1
          GO TO 999   !  Return
        ENDIF

!     REMEMBER CURRENT VERSION OF PP_STASH FOR NEXT TIME THRU LOOP
         OLD_STASH = PP_STASH(J)
!
      ENDDO

!     Check full length
      LEN = 0
      DO J=1,LEN2_LOOKUP
        LEN = LEN + MPP_DUMP_LEN(J)
      ENDDO

      IF (LEN  /=  LEN_DATA) THEN
        CMESSAGE = 'ADDR_CHK : Mismatch in length of data'
        WRITE (6,*) ' Length according to LOOKUP table ',LEN
        WRITE (6,*) ' Length set up in D1 array        ',LEN_DATA
        WRITE (6,*) ' You probably need to RECONFIGURE the start dump'
        ICODE = 2
        GO TO 999   !  Return
      ENDIF

 999  CONTINUE
      IF (lhook) CALL dr_hook('ADDRESS_CHECK',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ADDRESS_CHECK
