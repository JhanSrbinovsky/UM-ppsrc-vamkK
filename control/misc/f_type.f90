! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE F_TYPE-----------------------------------------------
!
!  Purpose:  Returns each field code and associated field length from
!            the PP header and a count of the number of fields
!            of each type.
!
! Programming standard : UMDP3
!
! Documentation: None
!
!-------------------------------------------------------------------
!
!*L  Arguments:-------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
      SUBROUTINE F_TYPE(LOOKUP,LEN2_LOOKUP,PP_NUM,N_TYPES               &
     &,PP_LEN,PP_STASH,PP_TYPE,PP_POS,PP_LS,FIXHD5,                      &
     &TITLE, verbose)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen
      USE Submodel_Mod
      IMPLICIT NONE

      INTEGER                                                           &
     & LEN2_LOOKUP                                                      &
                               !IN 2nd dimension of LOOKUP
     &,N_TYPES                                                          &
                               !IN No of separate field types in file
     &,LOOKUP(64,LEN2_LOOKUP)                                           &
                               !IN LOOKUP record
     &,PP_NUM(LEN2_LOOKUP)                                              &
                               !OUT No of successive fields with same co
     &,PP_LEN(LEN2_LOOKUP)                                              &
                               !OUT Length of field
     &,PP_STASH(LEN2_LOOKUP)                                            &
                               !OUT PP code of field
     &,PP_TYPE(LEN2_LOOKUP)                                             &
                               !OUT Integer/real/timeseries
     &,PP_POS(LEN2_LOOKUP)                                              &
                               !OUT Pointer to number of PP field
     &,PP_LS(LEN2_LOOKUP)                                               &
                               !OUT Data stored on land or sea pts
     &,FIXHD5                  !IN Fixed header item 5 (file type) 

      CHARACTER(LEN=80)TITLE
      LOGICAL :: verbose


! Local variables: -----------------------------------------------------
      INTEGER MODEL             !Internal model number from LOOKUP

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & PP_XREF(PPXREF_CODELEN)  !PPXREF codes for a given section/item

! External subroutines called:------------------------------------------
      CHARACTER(LEN=36) EXPPXC
      EXTERNAL EXPPXC
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------
      INTEGER                                                           &
     & ICODE                                                            &
                  ! Error code
     &,ITEM_CODE                                                        &
                  ! STASH item code
     &,SECTION    ! STASH section number

      CHARACTER                                                         &
     & CMESSAGE*80                                                      &
                   ! Error message
     &,PHRASE*(PPXREF_CHARLEN) ! Name of field

      INTEGER I,K
      character(len=20) :: valid(len2_lookup)
      integer :: fc_time(len2_lookup)
      integer :: lbc_levs(len2_lookup)
      integer :: lbproc(len2_lookup)
      integer :: lbtyp(len2_lookup)
      real,allocatable :: blev(:)
      real,allocatable :: bacc(:)
      
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!*----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('F_TYPE',zhook_in,zhook_handle)
      
      IF (verbose) THEN
        ALLOCATE(blev(len2_lookup))
        ALLOCATE(bacc(len2_lookup))
        blev(:)=0.0
        bacc(:)=0.0
      END IF

! Initialise arrays
      DO K=1,LEN2_LOOKUP
        PP_NUM(K)=1
        PP_LEN(K)=0
        PP_STASH(K)=0
        PP_TYPE(K)=0
        PP_POS(K)=0
        PP_LS(K)=0
        lbproc(k)=0
        lbtyp(k)=0
      ENDDO

      WRITE(valid(1),'(i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)')              &
     &                 lookup(1,1),':',lookup(2,1),':',                 &
     &                 lookup(3,1),':',lookup(4,1),':',                 &
     &                 lookup(5,1),':',lookup(6,1)
      fc_time(1)=lookup(14,1)

      N_TYPES=1
      PP_LEN(1)=LOOKUP(15,1)
      PP_STASH(1)=LOOKUP(42,1)
      lbproc(1)=lookup(25,1)
      lbtyp(1)=lookup(32,1)
      IF (verbose) THEN
        blev(1)=TRANSFER(lookup(52,1),blev(1))
        bacc(1)=TRANSFER(lookup(51,1),bacc(1))
      END IF
      PP_TYPE(1)=LOOKUP(39,1)
      PP_POS(1)=1
        IF(MOD(INT(LOOKUP(21,1)/10),10) == 2)THEN
          PP_LS(1)=MOD(INT(LOOKUP(21,1)/100),10)
        ENDIF

      DO K=2,LEN2_LOOKUP
        IF(LOOKUP(42,K) == LOOKUP(42,K-1).AND.                          &
     &     (LOOKUP(18,K)*lookup(19,k)) ==                               &
     &     (LOOKUP(18,K-1)*lookup(19,k-1))                              &
     &            .and.                                                 &
     &       lookup(1,k) == lookup(1,k-1) .and.                         &
     &       lookup(2,k) == lookup(2,k-1) .and.                         &
     &       lookup(3,k) == lookup(3,k-1) .and.                         &
     &       lookup(4,k) == lookup(4,k-1) .and.                         &
     &       lookup(5,k) == lookup(5,k-1) .and.                         &
     &       lookup(6,k) == lookup(6,k-1) .and.                         &
     &       lookup(25,k) == lookup(25,k-1) .and.                       &
     &       lookup(14,k) == lookup(14,k-1) .and.                       &
             .not. verbose )THEN
          PP_NUM(N_TYPES)=PP_NUM(N_TYPES)+1
        ELSE
          N_TYPES=N_TYPES+1
          PP_LEN(N_TYPES)=LOOKUP(15,K)
          PP_STASH(N_TYPES)=LOOKUP(42,K)
          lbproc(N_TYPES)=lookup(25,k)
          lbtyp(N_TYPES)=lookup(32,k)
          IF (verbose) THEN
            blev(N_TYPES)=TRANSFER(lookup(52,k),blev(N_TYPES))
            bacc(N_TYPES)=TRANSFER(lookup(51,k),bacc(N_TYPES))
          END IF
          PP_TYPE(N_TYPES)=LOOKUP(39,K)
          PP_POS(N_TYPES)=PP_POS(N_TYPES-1)+PP_NUM(N_TYPES-1)
          WRITE(valid(n_types),'(i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)')    &
     &lookup(1,k),':',lookup(2,k),':',                                  &
     &lookup(3,k),':',lookup(4,k),':',                                  &
     &lookup(5,k),':',lookup(6,k)
          fc_time(n_types)=lookup(14,k)
          if(fixhd5==5)then
            lbc_levs(n_types)=lookup(17,k)-100
          endif
        IF(MOD(INT(LOOKUP(21,K)/10),10) == 2)THEN
          PP_LS(N_TYPES)=MOD(INT(LOOKUP(21,K)/100),10)
        ENDIF
        ENDIF
      ENDDO

! Print out details of fields
      WRITE(6,'(/,A)') '************************************************'// &
                       '********************************'
      WRITE(6,'(''  '',/,'' '',A/)')TITLE

       SELECT CASE(fixhd5)
         CASE(3)
           WRITE(6,'(A)')'Key to output- fields file'
           WRITE(6,'(A)')'Data stored on land/sea points'
           WRITE(6,'(A)')'No of fields'
           IF (verbose) THEN
             WRITE(6,'(A)')'Level value'
             WRITE(6,'(A)')'Packing accuracy'
           END IF
           WRITE(6,'(A)')'Length of field before unpacking (==LBLREC)'
           WRITE(6,'(A)')'Data type'
           WRITE(6,'(A)')'STASH Code'
           WRITE(6,'(A)')'Met08 Code'
           WRITE(6,'(A)')'LBProc'
           WRITE(6,'(A)')'Start of first field'
           WRITE(6,'(A)')'Description'
           WRITE(6,'(A)')'Validity time (yyyy:mm:dd:hh:mn:ss)'
           WRITE(6,'(A)')'Forecast period'
         CASE(5)
           WRITE(6,'(A)')'Key to output - lbc file'
           WRITE(6,'(A)')'Data stored on land/sea points'
           WRITE(6,'(A)')'No of fields'
           WRITE(6,'(A)')'Length of field before unpacking (==LBLREC)'
           WRITE(6,'(A)')'Data type'
           WRITE(6,'(A)')'STASH Code'
           WRITE(6,'(A)')'Start of first field'
           WRITE(6,'(A)')'Description'
           WRITE(6,'(A)')'Validity time (yyyy:mm:dd:hh:mn:ss)'
         CASE default
           WRITE(6,'(A)')'Key to output - dump or other'
           WRITE(6,'(A)')'Data stored on land/sea points'
           WRITE(6,'(A)')'No of fields'
           WRITE(6,'(A)')'Length of field before unpacking (==LBLREC)'
           WRITE(6,'(A)')'Data type'
           WRITE(6,'(A)')'STASH Code'
           WRITE(6,'(A)')'Start of first field'
           WRITE(6,'(A)')'Description'
           WRITE(6,'(A)')'Validity time (yyyy:mm:dd:hh:mn:ss)'
           WRITE(6,'(A)')'Forecast period'
       END SELECT

      I=1
      DO K=1,N_TYPES
        if(lookup(42,i) == -99)then
          exit
        endif
        PHRASE=' '
        ITEM_CODE=MOD(LOOKUP(42,I),1000)
        SECTION=(LOOKUP(42,I)-ITEM_CODE)/1000
        MODEL=LOOKUP(45,I)
        ICODE = 0

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
        if(model == 10)then
          model = 1
        endif

! DEPENDS ON: exppxc
        PHRASE=EXPPXC(MODEL,SECTION,ITEM_CODE,                          &
     &              ICODE,CMESSAGE)
        IF(ICODE /= 0)THEN
          PHRASE='NON-STANDARD FIELD'
        ENDIF
        I=I+PP_NUM(K)
        SELECT CASE(fixhd5)
          CASE(3)
            IF(verbose) THEN
              WRITE(6,'('' '',I2,I5,F8.2,F6.1,I8,I2,4I6,1x,A36,1x,a22,i3)')    &
     &        PP_LS(K),PP_NUM(K),blev(k),bacc(k),PP_LEN(K),                  &
     &        PP_TYPE(K),PP_STASH(K),lbtyp(k),lbproc(K),PP_POS(K),           &
     &        PHRASE,valid(k),fc_time(k)
            ELSE
              WRITE(6,'('' '',I2,I5,I8,I4,4I6,1x,A36,1x,a22,i3)')         &
     &        PP_LS(K),PP_NUM(K),PP_LEN(K),                             &
     &        PP_TYPE(K),PP_STASH(K),lbtyp(k),lbproc(K),PP_POS(K),      &
     &        PHRASE,valid(k),fc_time(k)
            END IF
          CASE(5)
            WRITE(6,'('' '',I2,I5,I8,I4,2I6,1x,A36,1x,a22)')              &
     &      PP_LS(K),lbc_levs(K),PP_LEN(K),                             &
     &      PP_TYPE(K),PP_STASH(K),PP_POS(K),                           &
     &      PHRASE,valid(k)
          CASE default
            WRITE(6,'('' '',I2,I5,I8,I4,2I6,1x,A36,1x,a22,i3)')           &
     &      PP_LS(K),PP_NUM(K),PP_LEN(K),                               &
     &      PP_TYPE(K),PP_STASH(K),PP_POS(K),                           &
     &      PHRASE,valid(k),fc_time(k)
        END SELECT

      ENDDO
      
      IF (ALLOCATED(bacc)) DEALLOCATE(bacc)
      IF (ALLOCATED(blev)) DEALLOCATE(blev)
      WRITE(6,'(A,/)') '************************************************'// &
                       '********************************'

      IF (lhook) CALL dr_hook('F_TYPE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE F_TYPE
