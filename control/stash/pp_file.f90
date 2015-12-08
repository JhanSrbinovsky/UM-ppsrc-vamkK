! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE PP_FILE -----------------------------------------
!
!  Purpose:- To output a field to a PP_FILE
!
!  External documentation  C4
!
!--------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH     

SUBROUTINE pp_file(ppfield,icurrll,lenbuf,num_words,    &
    rmdi,comp_accrcy,pphoriz_out,unitpp,iwa,n_cols_out, &
    n_rows_out,packing,im_ident,                        &
    packing_type,current_io_pe,len_extra,               &
    srow_out,wcol_out,icode,cmessage)
  
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE stwork_aux, ONLY : &
      flushPendingStashForUnit
  USE IO, ONLY : &
      buffout,   &
      setpos
  USE IOS, ONLY :       &
      IOS_Setpos_Stash, &
      IOS_Write_Stash_PrePacked_Data
  USE IOS_STASH_COMMON, ONLY : &
      isUsingAsyncStash
  USE io_configuration_mod, ONLY : &
      io_field_padding

  USE UM_ParVars
  USE PrintStatus_mod,  ONLY : &
      printstatus,             &
      prstatus_diag

  USE Submodel_Mod

  USE chsunits_mod, ONLY : nunits

  IMPLICIT NONE
    
  CHARACTER(LEN=80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
  
  INTEGER, INTENT(IN) :: &
      LENBUF,       &! LENGTH OFF PP BUFFER
      UNITPP,       &! OUTPUT PP UNIT NUMBER
      im_ident,     &! Internal model identifier
      current_io_pe,&! PE which will do the I/O
      N_ROWS_OUT,   &! PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT
      N_COLS_OUT,   &! PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT
      COMP_ACCRCY,  &! PACKING ACCURACY IN POWER OF 2
      PPHORIZ_OUT,  &! SIZE OF OUTPUT FIELD
      ICURRLL,      &! Record Number
      LEN_EXTRA,    &! IN size of expected extra data
      SROW_OUT,     &! 1st southern row to output
      WCOL_OUT       ! 1st western column to output

  LOGICAL, INTENT(IN) :: &
      PACKING        !IN OVERALL Packing switch (T if pckng reqd)

  INTEGER, INTENT(OUT) :: &
      PACKING_TYPE   ! OUT set to 1 if WGDOS packing else set to zero.

  INTEGER, INTENT(INOUT) :: &
      NUM_WORDS,    &! NUMBER OF 64 BIT WORDS WORTH OF DATA
      ICODE          ! RETURN CODE FROM ROUTINE


  REAL                           &
      PPFIELD(PPHORIZ_OUT),      & !INOUT ARRAY TO STORE PPDATA
      BUFOUT(LENBUF+LEN_EXTRA),  & !OUTPUT PP BUFFER (ROUNDED UP)
      RMDI                         !IN     Missing data indicator


! Description: This include file contains information needed when
!              generating variable horizontal grid data in the
!              STASH extra data vector. Introduced UM 5.4 - R. Hill
!===================================================================
      LOGICAL :: X_VAR_GRID ! Whether variable grid in E-W direction
      LOGICAL :: Y_VAR_GRID ! and/or in S-N direction

      INTEGER :: VAR_GRID_TYPE ! 0 = none
                               ! 1 = T grid
                               ! 2 = U/V grid

      ! Grid boundaries for T and U,V
      REAL :: X_BOUNDARY(ROW_LENGTH_MAX+1,2)
      REAL :: Y_BOUNDARY(ROWS_MAX+1,2)

      ! Grid Points for T and U,V
      REAL :: X_GRID(ROW_LENGTH_MAX,2)
      REAL :: Y_GRID(ROWS_MAX,2)

      COMMON /OVARGRID/ X_VAR_GRID,Y_VAR_GRID                           &
     & ,X_BOUNDARY,Y_BOUNDARY,X_GRID,Y_GRID,VAR_GRID_TYPE

      ! The following parameters correspond to the extra data
      ! vector descriptors expected, for e.g., in PV-WAVE
      ! plotting routines (e.g. decode_extra.pro). There are
      ! numerous other areas of code where these integer
      ! descriptors must be handled (e.g. FIELDCOS, PPI2H, FTT)
      ! So it is not a trivial matter to introduce new code
      ! descriptors. Furthermore, ieee -32 will destroy these
      ! integers so PP data must always be processed via the
      ! long winded route: QXFIELDCOS -> PUTIBM -> FTT/PPI2H.
      ! (This thoroughly unsatisfactory state of affairs may
      ! be correctable with developments to ieee and convpp).
      INTEGER,PARAMETER :: x_coord_vector=1
                                     ! Indicates that an extra
                                     ! data vector gives LBNPT
                                     ! x-coordinate values
      INTEGER,PARAMETER :: y_coord_vector=2
                                     ! Indicates that an extra
                                     ! data vector gives LBROW
                                     ! Y-coordinate values
      INTEGER,PARAMETER :: x_lbnd_vector=12
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! x-boundary values
      INTEGER,PARAMETER :: x_ubnd_vector=13
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! x-boundary values
      INTEGER,PARAMETER :: y_lbnd_vector=14
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! y-boundary values
      INTEGER,PARAMETER :: y_ubnd_vector=15
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! y-boundary values
!dir$ cache_align bufout

!*  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
!*-------------------------------------------------------------------

!*-------------------------------------------------------------------
!   MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
!--------------------------------------------------------------------

!    DEFINE LOCAL VARIABLES
      INTEGER ::          &
          ml,             &! longitude counter
          jl,             &! latitude counter
          iwa,            &! record number
          ii,             &! counter
          ocode,          &! copy of the input value of icode
          length_fullwrd, &! length in bits of fullword var
          len_buf_words,  &! num_words rounded by 512 and actually
          num_out,        &! number of compressed (32 bit) words
          len_io,         &! not used, but needed for buffout call
          JJ               ! Local counter
      LOGICAL :: UV  
      REAL :: IX           ! RETURN VALUE FROM UNIT COMMAND

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('PP_FILE',zhook_in,zhook_handle)

!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) then
         goto 999
      ELSE IF (icode  <   0)then
         ocode = icode
         icode = 0
      END IF

!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
!======================================================================
      LENGTH_FULLWRD=64   !   LENGTH IN BITS OF FULLWORD VAR
!L   At this point packing,if required,will be done using the WGDOS
!L   method of packing.
      PACKING_TYPE=0
! Note the value of -26 corresponds to -15 (F) in ppxref.
! The packing acuracy is scaled to allow greater accuracy.
! Packing will only be attempted if there are at least 2 points per row
! in the PPfield.

      IF(PACKING .AND. N_COLS_OUT  >=  2) THEN
        ! 1. Climate wgdos packing has been selected for current file
        !    stream via UMUI
        IF(COMP_ACCRCY  >   -99) THEN
          ! 2. STASH packing profile for the field is set.
          PACKING_TYPE = 1
        ENDIF
      END IF
!
      IF(PACKING_TYPE == 1)THEN

! DEPENDS ON: coex
        CALL COEX(PPFIELD,PPHORIZ_OUT,BUFOUT,LENBUF,N_COLS_OUT,         &
     &  N_ROWS_OUT,NUM_OUT,COMP_ACCRCY,.TRUE.,RMDI,LENGTH_FULLWRD,      &
     &  icode,cmessage)

        NUM_WORDS=(NUM_OUT+1)/2 ! Round up to the nearest 64 Bit CRAY Wd
!  COEX returns the number of IBM words needed to hold the packed data
!                             ~~~
        LEN_BUF_WORDS=((NUM_WORDS+io_field_padding-1)/io_field_padding)*    &
     &   io_field_padding
IF (printstatus >= prstatus_diag) THEN
        WRITE(6,'(1X,A,I7,A,I7)')'PP_FILE: Field Packed from', pphoriz_out, &
        ' to', num_words
END IF
      ELSE IF(PACKING_TYPE == 4)THEN
        if (mype  ==  current_io_pe) then
! DEPENDS ON: runlen_encode
          CALL RUNLEN_ENCODE(PPFIELD,PPHORIZ_OUT,BUFOUT,PPHORIZ_OUT,    &
     &                     NUM_OUT,RMDI,ICODE,CMESSAGE)
        ! Size of run length encoded data is greater than unpacked
        ! field therefore leave field unpacked.
          if (NUM_OUT  >=  PPHORIZ_OUT) then
            PACKING_TYPE = 0
            DO JJ=1,PPHORIZ_OUT
              BUFOUT(JJ) = PPFIELD(JJ)
            END DO
            NUM_WORDS=PPHORIZ_OUT
            LEN_BUF_WORDS=LENBUF
          else
            NUM_WORDS=NUM_OUT
            LEN_BUF_WORDS=((NUM_WORDS+io_field_padding-1)/io_field_padding)*&
     &      io_field_padding
          endif

        end if
      ELSE  ! No packing required.
IF (printstatus >= prstatus_diag) THEN
        WRITE(6,'(1X,A,I7,I7)')'PP_FILE: Output field with grid dimensions', &
        N_ROWS_OUT,N_COLS_OUT 
END IF
        DO JJ=1,PPHORIZ_OUT
          BUFOUT(JJ) = PPFIELD(JJ)
        END DO

        NUM_WORDS=PPHORIZ_OUT
        LEN_BUF_WORDS=LENBUF
      ENDIF

      IF (VAR_GRID_TYPE >  0) THEN

         ! If we have a variable horizontal grid then we
         ! must add the grid data to the end of the PP record
         ! as "EXTRA DATA".
! DEPENDS ON: extra_variable_grid
         CALL EXTRA_VARIABLE_GRID(BUFOUT(NUM_WORDS+1),LEN_EXTRA         &
     &            ,SROW_OUT,N_ROWS_OUT                                  &
     &            ,WCOL_OUT,N_COLS_OUT)

         ! Adjust output buffer size for WFIO.
         NUM_WORDS = NUM_WORDS + LEN_EXTRA

         IF ((PACKING_TYPE == 1).OR.(PACKING_TYPE == 4)) THEN
           LEN_BUF_WORDS=((NUM_WORDS+io_field_padding-1)/io_field_padding)* &
     &                     io_field_padding ! No of words output
          ELSE
             LEN_BUF_WORDS=LENBUF+LEN_EXTRA
         ENDIF
      ENDIF ! If variable grid data needed adding

      IF (isUsingAsyncStash().AND.icurrll>=0) THEN
        !! If we are here we received a stash field object thats
        !  not using the async-queue method, and we are calling the 
        !  fallback code. Hence we should purge anything for the 
        !  unit that may be cached and not yet sent before we send 
        !  more data.
        CALL flushPendingStashForUnit(unitpp)
      END IF

      if (mype  ==  current_io_pe) then
      DO JJ=NUM_WORDS+1,LEN_BUF_WORDS
        BUFOUT(JJ)= 0.0
      ENDDO
      IF (isUsingAsyncStash().AND.icurrll>=0) THEN
        ! This is an alternate route for stash fields, where we use data
        ! on the backend IOS to keep a track of file locations. If we are 
        ! executing this code it means that some fields are written via async
        ! and this is a fallback for those that are not - in this case the 
        ! file locations are unknown and we need to use the stash aware setpos
        ! with the field number....
        CALL IOS_Setpos_Stash(unitpp,icurrll,icode)
        CALL IOS_Write_Stash_PrePacked_Data(unitpp,icurrll, &
            bufout(1:len_buf_words))
        len_io=len_buf_words
        ix = -1.0
      ELSE
        CALL setpos(unitpp,iwa,icode)
        CALL buffout(unitpp,bufout,len_buf_words,len_io,ix)
      END IF

!     WRITE(6,102) IWA,LEN_BUF_WORDS
  100 FORMAT(//,32X,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10F8.0/))
  102 FORMAT(' FROM PP_FILE    IWA  LEN_BUF_WORDS ',2I12)
!
      IF (IX /= -1.0.OR.LEN_IO /= LEN_BUF_WORDS) THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer out Data Field',IX,LEN_IO,                 &
     &                LEN_BUF_WORDS)
        CMESSAGE='PPFILE  : I/O error - PP Data Field Output'
        ICODE=7
        IF (lhook) CALL dr_hook('PP_FILE',zhook_out,zhook_handle)
        RETURN
      ENDIF
      endif ! (mype  ==  current_io_pe)
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF
  999 CONTINUE
      IF (lhook) CALL dr_hook('PP_FILE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PP_FILE
