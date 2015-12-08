! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read in a STASHmaster record

Module Rcf_ReadSTM_Mod

!  Subroutine Rcf_ReadSTM    - Reads in a STASHmaster record
!
! Description:
!   The routine reads in a STASHmaster record.
!
! Method:
!   A STASHmaster record (format as in UMDP C4) is read in on
!   unit nft and returned in STM_tmp.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   10/02/03   Extension of option code to 30 digits. T.White
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

SUBROUTINE Rcf_ReadSTM( STM_tmp, NFT, ICODE, CMESSAGE)

Use Rcf_Ppx_Info_mod, Only : &
    STM_record_type,     &
    STM_OptCodeLen,      &
    STM_VerCodeLen,      &
    Ppxref_pack_profs

IMPLICIT NONE

!     ARGUMENTS
INTEGER, INTENT(IN)             :: NFT        ! UNIT NUMBER STMSTS

INTEGER, INTENT(OUT)            :: ICODE      ! RETURN CODE
CHARACTER (Len=80), Intent(Out) :: CMESSAGE   ! RETURN MESSAGE IF THERE
                                              ! IS A FAILURE
Type (STM_record_type), Intent(Out) :: STM_tmp ! returned data

! Local Variables
INTEGER       :: IMSK          ! Decimal equivalent of binary IMASK
INTEGER       :: II            ! loop mark.
Integer       :: jj
Integer       :: opcod( STM_OptCodeLen )     ! option codes
Integer       :: imask( STM_VerCodeLen )     ! version mask
! Format spec. for variable length codes
Character(Len=27) :: CodeFormat

write( CodeFormat, '(A,I2,A,I2,A)' ) &
       '(3X,',                       &
       STM_OptCodeLen,               &
       '(I1),3X,',                   &
       STM_VerCodeLen,               &
       '(I1),2X,I5)'

ICODE=0
CMESSAGE=' '

READ(NFT,2010,END=3100,ERR=3200) &
      STM_tmp % model,           &
      STM_tmp % section,         &
      STM_tmp % item,            &
      STM_tmp % name
2010 FORMAT(2X,3(I5,2X),A36)

IF (STM_tmp % model == -1) GO TO 9999

READ(NFT,2110,END=3100,ERR=3200) &
      STM_tmp % space_code,      &
      STM_tmp % ptr_code,        &
      STM_tmp % timavail_code,   &
      STM_tmp % grid_type,       &
      STM_tmp % lv_code,         &
      STM_tmp % lb_code,         &
      STM_tmp % lt_code,         &
      STM_tmp % pt_code,         &
      STM_tmp % pf_code,         &
      STM_tmp % pl_code,         &
      STM_tmp % lev_flag
2110 FORMAT(2X,11(I5,2X))

Read( NFT, CodeFormat, End=3100, Err=3200 ) &
      (opcod(ii) , ii=1,STM_OptCodeLen),    &
      (Imask(ii) , ii=1,STM_VerCodeLen),    &
      STM_tmp % halo_type

! Option code converted to 5-digit groups
Do ii = 1, STM_OptCodeLen / 5
   jj = STM_OptCodeLen - (ii-1)*5
 ! ii = 1, 2, 3, 4, ...
 ! jj = 30, 25, 20, ...
   STM_tmp % opt_code(ii) =                             &
        opcod(jj) + opcod(jj-1)*10   + opcod(jj-2)*100  &
                  + opcod(jj-3)*1000 + opcod(jj-4)*10000
EndDo

!   Binary version mask was read into array IMASK
!   Convert version mask to decimal form IMSK
    IMSK = 0
    DO II=20,1,-1
      IF((IMASK(II).NE.0).AND.(IMASK(II).NE.1)) THEN
        WRITE(6,*) 'Rcf_ReadSTM: improper IMASK in user diag'
        WRITE(6,*) 'Model, Section, Item ', &
             STM_tmp % model,               &
             STM_tmp % section,             &
             STM_tmp % item
      ELSE
        IF(IMASK(II).EQ.1) THEN
          IMSK=IMSK+2**(20-II)
        END IF
      END IF
    END DO
!     Insert decimal value of version mask
    STM_tmp % version_mask = IMSK

READ(NFT,2130,END=3100,ERR=3200) &
      STM_tmp % data_type,       &
      STM_tmp % dump_packing,    &
     (STM_tmp % packing_acc(II), &
      II= 1, PPXREF_PACK_PROFS)
2130 FORMAT(2X,I5,2X,I5,3X,I3,9(2X,I3))

READ(NFT,2140,END=3100,ERR=3200) &
      STM_tmp % rotate_code,     &
      STM_tmp % field_code,      &
      STM_tmp % user_code,       &
      STM_tmp % lbvc_code,       &
      STM_tmp % base_level,      &
      STM_tmp % top_level,       &
      STM_tmp % ref_lbvc_code,   &
      STM_tmp % cf_levelcode,    &
      STM_tmp % cf_fieldcode
2140 FORMAT(2X,9(I5,2X))

3100 GO TO 9999 ! Normal completion
3200 WRITE(6,*)' MESSAGE FROM ROUTINE Rcf_ReadSTM: '
WRITE(6,*)' ERROR OCCURRED WHILE READING STASHmaster FILE '
CMESSAGE=' Rcf_ReadSTM: ERROR READING STASHMASTERS FILE'
ICODE=2

9999 CONTINUE
RETURN

END Subroutine Rcf_ReadSTM

End Module Rcf_ReadSTM_Mod
