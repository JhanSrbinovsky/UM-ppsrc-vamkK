! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Get Ancillary file records from ANCILfiles

Module GetAnc_Files_Mod

!  Subroutine GetAnc_Files  - Get Ancillary file records
!
! Description:
!    Get records from ANCILfiles files and insert into table.
!
! Method:
!    1. Opens the ANCILfiles file and reads in single records until
!       an end of file record is read in.
!    2. For system files, the records are added to the table
!       of records.
!    3. For User files, the record will overwrite an existing
!       record otherwise it is added to the table of records.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  D. Robinson
!   5.3   11/10/01   Alter to use Max_Filename_Len parameter. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine GetAnc_Files ( AncMsrNam, EnvDir, UsrANC )

Use Submodel_Mod, Only :      &
    Internal_Model_Index

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    Rcf_Free_Unit,            &
    Max_Filename_Len

Use Ereport_Mod, Only :    &
    Ereport

Use Ancil_Mod, Only :      &
    ancFiles,              &
    max_ancFiles,          &
    anc_file_type,         &
    anc_file

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

Use ReadAnc_Mod, Only :   &
    ReadAnc_Files

Use BcastANC_Mod, Only :  &
    BcastANC_Files

USE UM_ParVars, Only :   &
    mype,                 &
    nproc

IMPLICIT NONE

! Arguments
Character (Len=13):: AncMsrNam   ! Names of ancil master file
Character (Len=9) :: EnvDir      ! Env variable for directory name
Logical, Optional :: UsrANC      ! Logical controlling User Ancilmaster

! Local Variables
Integer            :: I             ! Loop counters
Integer            :: Irec          !
Integer            :: ErrorStatus   ! Error return code
Integer            :: nftancilm     ! Unit number for ANCILmaster
Integer            :: IOSTATUS
Integer            :: LAncFileNo
Integer            :: UAncRow
Integer            :: FirstBlank    ! Used to append filename to dir
Integer            :: Info          ! For GCOM
Integer            :: msg           ! For GCOM tag
Integer            :: File_Len
Integer            :: Anc_File_No

Character (Len=Max_Filename_Len) :: ANCIL_MSTR
                                    ! Full pathname for ANCIL
                                    ! master files
Character (Len=80) :: CMESSAGE      ! Error return message
Character (Len=1)  :: CHAR1
Character (Len=*), Parameter :: RoutineName='Getanc_Files'

Logical            :: OVERWRITE     ! T if an ancilmaster record is
                                    ! being overwritten by a user rec
Logical            :: UUsrANC       ! Controlling UserSTASH

Type (ANC_file_type) :: ANC_tmp   ! Temporary record

! --------------------------------------------------------------

ErrorStatus = 0
IOStatus    = 0

! Set default for UserANCIL flag
If ( Present( UsrANC ) ) Then
  UUsrANC = UsrANC
Else
  UUsrANC = .False.
Endif

! Allocate space for the ANCILmaster records.
! If one needs allocating, they all do.
If ( .Not. Allocated ( anc_file ) ) Then
  Allocate( anc_file ( max_ancFiles ) )

! Initialise all ancillary file numbers to zero
  ANC_file(:) % anc_file_number = 0

End If

!--------------------------------
!Read in records from ANCILmaster
!--------------------------------

If ( mype == 0 ) Then

  ! Get a file unit
  Call Rcf_Get_Unit( nftancilm )

  ! Get directory name for ANCILmaster
  Call Fort_Get_Env( EnvDir, 9, ANCIL_MSTR, Max_Filename_Len,      &
                     ErrorStatus )

  If ( ErrorStatus /= 0 ) Then
    Cmessage = 'Cannot get ANCIL directory from Env. Vars'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  Endif

! Append filename
  FirstBlank = Len_Trim( ANCIL_MSTR ) + 1
  ANCIL_MSTR = Trim( ANCIL_MSTR )

  File_Len = Len_Trim (AncMsrNam)

  ANCIL_MSTR(FirstBlank:FirstBlank)='/'
  ANCIL_MSTR(FirstBlank+1:FirstBlank+File_Len)=AncMsrNam

! Existence of file already checked in hdancrf

! Open the ANCILmaster file
  Open( Unit=NFTANCILM, File=ANCIL_MSTR, Iostat=IOStatus )

  If (IOStatus.NE.0) Then
    Write (6,*) 'ERROR in routine GETANC_FILES'
    Write (6,*) &
    'CANNOT OPEN ANCILmaster FILE, IOSTATUS=',IOStatus
    Write (6,*) 'UNIT=',Nftancilm,' FILE=',ANCIL_MSTR
    ErrorStatus=100
    CMESSAGE='ERROR OPENING ANCILmaster'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If   ! pe0

Anc_File_No = 0
Do While (Anc_File_No /= -1)

  If ( mype == 0 ) Then
    Read( nftancilm, '(A1)' ) CHAR1
  End If

  msg = 9051
  Call gc_cbcast( msg, 1, 0, nproc, info, char1 )

  If (CHAR1.EQ.'1') Then

    !Read block of records
    If ( mype == 0 ) Then
      Backspace NFTANCILM
      Call ReadAnc_Files( ANC_tmp, NFTANCILM, ErrorStatus, CMESSAGE)
    End If

    Call BcastANC_Files ( ANC_tmp, 0 )

    Anc_file_No = ANC_tmp % anc_file_number

    If (PrintStatus > PrStatus_Diag .and. mype == 0 ) then
      If (Anc_File_No /= -1) then

      write (6,*) ' ====================================='
      write (6,*) ' ancil_fn ',ANC_tmp % anc_file_number
      write (6,*) ' model_no ',ANC_tmp % model_number
      write (6,*) ' ancil_ev ',ANC_tmp % anc_env_var
      write (6,*) ' ancil_ft ',ANC_tmp % anc_file_title
      write (6,*) ' ====================================='

      endif
    endif

    If (Anc_File_No /= -1) Then      ! Not end of file

      If (.not. UUsrANC) Then  !  Not an User ANCILmaster file

        !   Increment row number
        ancFiles = ancFiles + 1   !   Increment row number

        !   Insert new record
        anc_file( ancFiles ) = ANC_tmp

        !   Set flag to indicate record came from ANCILmaster file
        anc_file( ancFiles ) % anc_flag = 'A'

        If (ancFiles .GT. max_ancFiles) Then
          ErrorStatus = 10
          Write ( Cmessage , '(A, I4)' ) &
          ' Too many Ancil records read in ',ancFiles
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

      Else      !  User ANCILmaster record

        !   Initialise OVERWRITE  switch
        OVERWRITE  = .FALSE.
        char1     = '0'

        !   Insert user record.
        !   No. of records extracted from ANCILmaster file(s)= ancFiles
        UAncRow    =   0

        Do I = 1, max_ancFiles

          !  Get ancil_file _no for this row
          LAncFileNo = anc_file ( I ) % anc_file_number
!         LModel     = 1

          If ( PrintStatus > PrStatus_Diag .and. mype==0 ) then
            write (6,*) 'Row No (I) ',I
            write (6,*) ' LAncFileNo   ',LAncFileNo
            write (6,*) ' Anc_File_no  ',anc_file_no
            write (6,*) ' UAncRow  ',UAncRow
          endif

          If (LAncFileNo > 0 .AND. UAncRow == 0) THEN

            ! Check whether previous item is being overwritten
            If (LAncFileNo == Anc_File_No ) Then

              If (anc_file(I) % anc_flag == 'A') Then

                OVERWRITE =.TRUE.

                If ( PrintStatus >= PrStatus_Normal .and. &
                     mype == 0 ) Then
                  Write (6,*) &
                  ' ANCILmaster record ',I,' has been overwritten by'
                  Write (6,*) &
                  ' User ANCILmaster record '
                  Write (6,*) ' Anc Ref No ',Anc_File_No
                End If

              Else If (anc_file(I) % anc_flag == 'U') Then

                Write (6,*) ' Error in GETANC_FILES.'
                Write (6,*) ' User diagnostic duplicated'
                Write (6,*) ' Anc File No ',Anc_File_No
                ErrorStatus=100
                CMESSAGE=' User diagnostic duplicated'
                Call Ereport( RoutineName, ErrorStatus, Cmessage )

              End If

            End If

            !     Determine appropriate row number
            If (LAncFileNo == Anc_File_No .AND. &
                UAncRow   == 0)             THEN

              UAncRow=I    ! Row number found

              !  Over-write pre-existing record
              anc_file( UAncRow ) = ANC_tmp

              !  Set flag to indicate User record
              anc_file( UAncRow ) % anc_flag = 'U'

            Else If (LAncFileNo > Anc_File_No .and. UAncRow == 0 ) Then

              UAncRow=I    ! Row number found

              ! This record will be inserted between two pre-existing
              ! records. Create a spare row by moving all subsequent
              ! records down.

              write (6,*) ' User ANCILmaster record to be '
              write (6,*) ' inserted between two pre-existing records.'
              write (6,*) ' Ancil File No ',Anc_File_No

              write (6,*) ' Moving records ',UAncRow,':',AncFiles,&
                        ' to ',UAncRow+1,':',AncFiles+1

              ! Move records down to create spare row
              Do irec = ancFiles+1, UAncRow+1, -1
                anc_file( irec ) = anc_file( irec - 1 )
              End Do

              !  Insert new record
              anc_file( UAncRow ) = ANC_tmp

              !  Set flag to indicate User record
              anc_file( UAncRow ) % anc_flag = 'U'

            End If

          Else If (Anc_file_No == 0 .AND. UAncRow == 0) Then

            ! This record will be added after all pre-existing records

            write (6,*) ' User ANCILmaster record to be '
            write (6,*) ' inserted after all pre-existing records.'
            Write (6,*) ' Anc File No ',Anc_File_No

            UAncRow = I

            !  Add new record
            anc_file( UAncRow ) = ANC_tmp

            !  Set flag to indicate User record
            anc_file( UAncRow ) % anc_flag = 'U'

          End If

        End Do   !  Loop over i

        ! Increment ancFiles if UserANCIL record has been added.
        If (.NOT. OVERWRITE ) Then

          ancFiles = ancFiles + 1

!         Remove eventually
          If (PrintStatus >= PrStatus_Diag .and. mype == 0) then
            write (6,*) ' ancFiles increased to ',ancFiles
            write (6,*) ' Ancil File Nos : ', &
            anc_file(:) % anc_file_number
          endif

        End If  !  If .not.overwrite

      End If  !  If UAncFileNo > 0

    End If  !  If .not.UsrAnc

  End If  !  If Anc_File_No /= -1

End Do  !   Loop over records until Anc_File No = -1 raeched

! Close relevant file and get rid of unit
If ( mype == 0 ) Then
  Close( nftancilm )
  Call Rcf_Free_Unit( nftancilm )
End If

Return
End Subroutine GetAnc_Files

End Module GetAnc_Files_Mod
