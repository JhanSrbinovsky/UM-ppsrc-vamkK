! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Module to read ancil-master file records

Module ReadAnc_Mod

!  Subroutine ReadAnc_Fields      - reads ancil-master FIELD record
!  Subroutine ReadAnc_Files       - reads ancil-master FILE  record
!
! Description:
!    The module reads in an ancil-master file record
!
! Method:
!    Data is read in on Unit NFT and returned in the ANC_TMP record.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  D.Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

!------------------------------------------------------------------
! Routine to read in ancil-master record for ancillary fields
!------------------------------------------------------------------

SUBROUTINE ReadAnc_Fields ( ANC_tmp, NFT, ICODE, CMESSAGE)

Use ancil_mod, Only : &
    anc_record_type

Implicit None

!     Arguments
Type (ANC_record_type), Intent(Out)  :: ANC_tmp  ! returned data
INTEGER, INTENT(IN)             :: NFT    ! UNIT NUMBER STMSTS
INTEGER, INTENT(OUT)            :: ICODE  ! RETURN CODE
CHARACTER (Len=80), Intent(Out) :: Cmessage  ! RETURN MESSAGE IF
                                             ! THERE IS A FAILURE

ICODE=0
CMESSAGE=' '

READ(NFT,2010,END=3100,ERR=3200) &
      ANC_tmp % ancil_ref_number  , &
      ANC_tmp % model_number      , &
      ANC_tmp % section_number    , &
      ANC_tmp % item_number       , &
      ANC_tmp % anc_name
2010 FORMAT(2X,4(I5,2X),A36)

IF (ANC_tmp % model_number == -1) GO TO 9999

!READ(NFT,2110,END=3100,ERR=3200) &
!      ANC_tmp % anc_file_number,   &
!      ANC_tmp % anc_env_var,        &
!      ANC_tmp % anc_file_title
! 2110 FORMAT(2X,I5,2X,A8,1X,A41)

READ(NFT,2110,END=3100,ERR=3200) &
      ANC_tmp % anc_file_number
2110 FORMAT(2X,I5)


3100 GO TO 9999 ! Normal completion
3200 WRITE(6,*)' MESSAGE FROM ROUTINE ReadAnc: '
WRITE(6,*)' ERROR OCCURRED WHILE READING ANCILmaster FILE '
CMESSAGE=' ReadAnc: ERROR READING ANCILmaster FILE'
ICODE=2

9999 CONTINUE
RETURN

END Subroutine ReadAnc_Fields

!------------------------------------------------------------------
! Routine to read in ancil-master record for ancillary files
!------------------------------------------------------------------

SUBROUTINE ReadAnc_Files ( ANC_tmp, NFT, ICODE, CMESSAGE)

Use ancil_mod, Only : &
    anc_file_type

Implicit None

! Arguments
Type (anc_file_type), Intent(Out)  :: ANC_tmp  ! returned data
Integer, Intent(IN)             :: NFT    ! UNIT NUMBER STMSTS
Integer, Intent(OUT)            :: ICODE  ! RETURN CODE
Character (Len=80), Intent(Out) :: Cmessage  ! RETURN MESSAGE IF
                                             ! THERE IS A FAILURE

ICODE=0
CMESSAGE=' '

Read (Nft,2010, END=3100, ERR=3200)  &
      ANC_tmp % anc_file_number    , &
      ANC_tmp % model_number       , &
      ANC_tmp % anc_env_var        , &
      ANC_tmp % anc_file_title
! 2010 FORMAT(2X,I5,2X,A8,2X,A41)
2010 FORMAT(2X,I5,2X,I4,2X,A8,1X,A41)

IF (ANC_tmp % anc_file_number == -1) GO TO 9999

3100 GO TO 9999 ! Normal completion
3200 WRITE(6,*)' MESSAGE FROM ROUTINE ReadAnc_Files: '
WRITE(6,*)' ERROR OCCURRED WHILE READING ANCILmaster FILE '
CMESSAGE=' ReadAnc: ERROR READING ANCILmaster FILE'
ICODE=2

9999 CONTINUE
Return
End Subroutine ReadAnc_Files


End Module ReadAnc_Mod
