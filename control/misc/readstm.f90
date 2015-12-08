! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE READSTM
!   PURPOSE  TO READ A RECORD FROM THE PRE STASH-MASTERS FILE
! AND RETURN THE PPXREF CODES AND NAME OF A GIVEN DIAGNOSTIC
! 
! LOGICAL COMPONENT R913
! 
! PROJECT TASK: C4
! 
! PROGRAMMING STANDARD  UMDP 4
! 
! EXTERNAL DOCUMENT C4
!
!
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Misc
      SUBROUTINE READSTM                                                &
     &          (IMASK,DNAM,CODES,NFT,ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE cppxref_mod
      USE Submodel_Mod
      IMPLICIT NONE


!     ARGUMENTS
      INTEGER NFT             !IN:    UNIT NUMBER STMSTS

      INTEGER ICODE           !OUT: RETURN CODE
      CHARACTER(LEN=80) CMESSAGE !OUT: RETURN MESSAGE IF THERE IS A FAILURE

      CHARACTER(LEN=1) DNAM(PPXREF_CHARLEN) !OUT: VARIABLE NAME FROM RECORD
      INTEGER CODES(PPXREF_CODELEN)    !OUT: PPXREF CODES FROM RECORD
      INTEGER IMASK(20)                !OUT: VERSION MASK

      INTEGER IMSK          ! Decimal equivalent of binary IMASK
      INTEGER II                       !LOCAL: loop mark.
      integer opcod(30)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
      IF (lhook) CALL dr_hook('READSTM',zhook_in,zhook_handle)
      ICODE=0
      CMESSAGE=' '

      READ(NFT,2010,END=3100,ERR=3200)                                  &
     & CODES(ppx_model_number)  ,                                       &
     & CODES(ppx_section_number),                                       &
     & CODES(ppx_item_number)   , DNAM
 2010 FORMAT(2X,3(I5,2X),36A1)

      IF(CODES(ppx_model_number) == -1) GO TO 9999

      READ(NFT,2110,END=3100,ERR=3200)                                  &
     & CODES(ppx_space_code),                                           &
     & CODES(ppx_ptr_code),                                             &
     & CODES(ppx_timavail_code),                                        &
     & CODES(ppx_grid_type),                                            &
     & CODES(ppx_lv_code),                                              &
     & CODES(ppx_lb_code),                                              &
     & CODES(ppx_lt_code),                                              &
     & CODES(ppx_pt_code),                                              &
     & CODES(ppx_pf_code),                                              &
     & CODES(ppx_pl_code),                                              &
     & CODES(ppx_lev_flag)
 2110 FORMAT(2X,11(I5,2X))

! 30-digit option code read in as 6x5 digit groups instead of 4x5
      READ(NFT,2120,END=3100,ERR=3200)                                  &
     &(opcod(ii),ii=1,30),                                              &
     &(IMASK(II),II=1,20),                                              &
     &CODES(ppx_halo_type)
 2120 FORMAT(3X,30(I1),3X,20(I1),2X,I5)

      CODES(ppx_opt_code  )=                                            &
     & opcod(30)+opcod(29)*10  +opcod(28)*100  +                        &
     &           opcod(27)*1000+opcod(26)*10000
      CODES(ppx_opt_code+1)=                                            &
     & opcod(25)+opcod(24)*10  +opcod(23)*100  +                        &
     &           opcod(22)*1000+opcod(21)*10000
      CODES(ppx_opt_code+2)=                                            &
     & opcod(20)+opcod(19)*10  +opcod(18)*100  +                        &
     &           opcod(17)*1000+opcod(16)*10000
      CODES(ppx_opt_code+3)=                                            &
     & opcod(15)+opcod(14)*10  +opcod(13)*100  +                        &
     &           opcod(12)*1000+opcod(11)*10000
      CODES(ppx_opt_code+4)=                                            &
     & opcod(10)+opcod( 9)*10  +opcod( 8)*100  +                        &
     &           opcod( 7)*1000+opcod( 6)*10000
      CODES(ppx_opt_code+5)=                                            &
     & opcod( 5)+opcod( 4)*10  +opcod( 3)*100  +                        &
     &           opcod( 2)*1000+opcod( 1)*10000

!   Binary version mask was read into array IMASK
!   Convert version mask to decimal form IMSK
          IMSK = 0
          DO II=20,1,-1
            IF((IMASK(II) /= 0).AND.(IMASK(II) /= 1)) THEN
              WRITE(6,*) 'READSTM: improper IMASK in user diag'
              WRITE(6,*) 'Model, Section, Item ',                       &
     &        CODES(ppx_model_number)  ,                                &
     &        CODES(ppx_section_number),                                &
     &        CODES(ppx_item_number)
            ELSE
              IF(IMASK(II) == 1) THEN
                IMSK=IMSK+2**(20-II)
              END IF
            END IF
          END DO
!     Insert decimal value of version mask
          CODES(ppx_version_mask)=IMSK

      READ(NFT,2130,END=3100,ERR=3200)                                  &
     & CODES(ppx_data_type),                                            &
     & CODES(ppx_dump_packing),                                         &
     &(CODES(II),                                                       &
     & II= ppx_pack_acc, ppx_pack_acc+PPXREF_PACK_PROFS-1)
 2130 FORMAT(2X,I5,2X,I5,3X,I3,9(2X,I3))

      READ(NFT,2140,END=3100,ERR=3200)                                  &
     & CODES(ppx_rotate_code),                                          &
     & CODES(ppx_field_code),                                           &
     & CODES(ppx_user_code),                                            &
     & CODES(ppx_lbvc_code),                                            &
     & CODES(ppx_base_level),                                           &
     & CODES(ppx_top_level),                                            &
     & CODES(ppx_ref_lbvc_code),                                        &
     & CODES(ppx_cf_levelcode),                                         &
     & CODES(ppx_cf_fieldcode)
 2140 FORMAT(2X,9(I5,2X))
 3100 GO TO 9999 ! Normal completion
 3200 WRITE(6,*)' MESSAGE FROM ROUTINE READSTM: '
      WRITE(6,*)' ERROR OCCURRED WHILE READING STASHmaster FILE '
      CMESSAGE=' READSTM: ERROR READING STASHMASTERS FILE'
      ICODE=2

 9999 CONTINUE
      IF (lhook) CALL dr_hook('READSTM',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE READSTM

