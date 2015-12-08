! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Get Ancillary Field records from ANCILfields

Module GetAnc_Fields_Mod

!  Subroutine GetAnc_Fields  - Get Ancillary field records
!
! Description:
!    Get records from ANCILfields files and insert into table.
!
! Method:
!    1. Opens the ANCILfields file and reads in single records until
!       an end of file record is read in.
!    2. For system files, the records are added to the table
!       of records.
!    3. For User files, the record will overwrite an existing
!       record otherwise it is added to the table of records.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine GetAnc_Fields ( AncMsrNam, EnvDir, UsrANC )

Use Submodel_Mod, Only :      &
    Internal_Model_Index

Use Rcf_FortranIO_Mod, Only :     &
    Rcf_Get_Unit,             &
    Rcf_Free_Unit,            &
    Max_Filename_Len

Use Ereport_Mod, Only :   &
    Ereport

Use Ancil_Mod, Only :      &
    ancrecs,               &
    max_ancrecs,           &
    ANC_record_type,       &
    ANC_record

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

Use ReadAnc_Mod, Only :   &
    ReadAnc_Fields

Use BcastANC_Mod, Only :  &
    BcastANC_Fields

USE UM_ParVars, Only :   &
    mype,                 &
    nproc

IMPLICIT NONE

! Arguments
Character (Len=13):: AncMsrNam   ! Names of ancil master files
Character (Len=9) :: EnvDir      ! Env variable for directory name
Logical, Optional :: UsrANC      ! Logical controlling User Ancilmaster

! Local Variables
Integer            :: I, ID         ! Loop counters
Integer            :: ErrorStatus   ! Error return code
Integer            :: nftancilm     ! Unit number for ANCILmaster file
Integer            :: IOSTATUS
Integer            :: Im_ident      !
Integer            :: Anc_Ref_No    !
Integer            :: Section       !
Integer            :: Item          !
Integer            :: LModel
Integer            :: LAncRefNo
Integer            :: UAncRow
Integer            :: FirstBlank    !Used to append Upsm filename to dir
Integer            :: Info          ! For GCOM
Integer            :: msg           ! For GCOM tag
Integer            :: File_Len

Character (Len=Max_Filename_Len) :: ANCIL_MSTR
                                    ! Full pathname for ANCIL
                                    ! master files
Character (Len=80) :: CMESSAGE      ! Error return message
Character (Len=1)  :: CHAR1
Character (Len=*), Parameter :: RoutineName='Getanc_Fields'

Logical            :: OVERWRITE     ! T if an ancil master record is
                                    ! being overwritten by a user rec
Logical            :: UUsrANC       ! Controlling User AncilMaster

Type (ANC_record_type) :: ANC_tmp   ! Temporary record

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
If ( .Not. Allocated ( ANC_record ) ) Then
  Allocate( ANC_record( max_ancRecs ) )

! Initialise
  ANC_record(:) % ancil_ref_number = 0
  ANC_record(:) % model_number = 0

End If

!--------------------------------
!Read in records from ANCILmaster
!--------------------------------
! Get a file unit

If ( mype == 0 ) Then

! Get a file unit
  Call Rcf_Get_Unit( nftancilm )

! Get directory name for ANCILmaster
  Call Fort_Get_Env( EnvDir, 9, ANCIL_MSTR, Max_Filename_Len,     &
                     ErrorStatus )

  If ( ErrorStatus /= 0 ) Then
    Cmessage = 'Cannot get ANCIL directory from Env. Vars'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  Endif

! Append filename
  firstblank = LEN_TRIM( ancil_mstr ) + 1
  ancil_mstr = TRIM( ancil_mstr )

  file_len = LEN_TRIM(ancmsrnam)
  IF (uusranc) THEN
    ancil_mstr(firstblank:firstblank)='.'
  ELSE
    ancil_mstr(firstblank:firstblank)='/'
  END IF

  ancil_mstr(firstblank+1:firstblank+file_len)=ancmsrnam

! Open the ANCILmaster file
  Open( Unit=nftancilm, File=ANCIL_MSTR, Iostat=IOStatus )

  If (IOStatus /= 0) Then
    Write (6,*) 'ERROR in routine GETANC_FIELDS'
    Write (6,*) &
    ' CANNOT OPEN ANCILmaster FILE, IOSTATUS=',IOStatus
    Write (6,*) 'UNIT=',NFTANCILM,' FILE=',ANCIL_MSTR
    ErrorStatus=100
    CMESSAGE=' GETANC: ERROR OPENING ANCILmaster'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If   ! pe0

Im_ident = 0
Do While (Im_ident /= -1)

  If ( mype == 0 ) Then
    Read( nftancilm, '(A1)' ) CHAR1
  End If

  msg = 9051
  Call gc_cbcast( msg, 1, 0, nproc, info, char1 )

  If (CHAR1.EQ.'1') Then

    !Read block of records
    If ( mype == 0 ) Then
      Backspace NFTANCILM
      Call ReadAnc_Fields( ANC_tmp, NFTANCILM, ErrorStatus, CMESSAGE)
    End If

    Call BcastANC_Fields ( ANC_tmp, 0 )

    Anc_Ref_No = ANC_tmp % ancil_ref_number
    Im_ident   = ANC_tmp % model_number
    Item       = ANC_tmp % item_number
    Section    = ANC_tmp % section_number
    
    If (PrintStatus > PrStatus_Diag .and. mype == 0 ) then
      If (Im_ident /= -1) then

        write (6,*) ' ====================================='
        write (6,*) ' anc_ref_no ',anc_ref_no
        write (6,*) ' model    ',im_ident
        write (6,*) ' section  ',section
        write (6,*) ' stashc   ',item
        write (6,*) ' anc_name ',ANC_tmp % anc_name
        write (6,*) ' ancil_fn ',ANC_tmp % anc_file_number
        write (6,*) ' ====================================='

      endif
    endif

    If (Im_ident /= -1) Then      ! Not end of file

      If (.not. UUsrANC) then  !  Not an User ANCILmaster file

        !  Increment row number
        ancRecs = ancRecs + 1

        !  Insert new record
        ANC_record( ancRecs ) = ANC_tmp

        !  Set flag to indicate record came from ANCILmaster file
        ANC_record( ancRecs ) % anc_flag = 'A'

        If (ancRecs .GT. max_ancRecs) Then
          ErrorStatus = 10
          Write ( Cmessage , '(A, I4)' ) &
          ' Too many Ancil records read in: ',ancRecs
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

      Else   !  User ANCILmaster file

        !   Initialise OVERWRITE  switch
        OVERWRITE  = .FALSE.

        !   Insert user record.
        !   No. of records extracted from ANCILmaster file(s)= ancRecs.
        UAncRow    =   0

        Do I = 1, max_ancRecs

          !     Determine values of model and ancil_ref_no for this row
          LAncRefNo =   ANC_record( I ) % ancil_ref_number
          LModel    =   ANC_record( I ) % model_number

          If ( PrintStatus > PrStatus_Diag .and. mype==0 ) then
            write (6,*) 'Row No (I) ',I
            write (6,*) ' LAncRefNo   ',LAncRefNo, ' LModel ',LModel
            write (6,*) ' Anc_ref_no  ',anc_ref_no,' Model  ',im_ident
            write (6,*) ' UAncRow  ',UAncRow
          endif

            If (LAncRefNo > 0 .AND. LModel > 0 .AND. UAncRow == 0) THEN

              !     Check whether previous item is being overwritten
              If (Im_ident == LModel   .AND. &
              Anc_Ref_No == LAncRefNo ) THEN

                If      (ANC_record(I) % anc_flag == 'A') Then

                  OVERWRITE =.TRUE.
                  If ( mype == 0 .and. &
                       PrintStatus >= PrStatus_Normal) Then
                    Write (6,*) &
                    ' ANCILmaster record ',I,' has been overwritten by'
                    Write (6,*) &
                    ' User ANCILmaster record '
                    Write (6,*) ' Anc Ref No ',Anc_Ref_No, &
                    ' Model ',Im_ident
                  End If

                Else If (ANC_record(I) % anc_flag == 'U') Then

                  Write (6,*) ' Error in GETANC_FIELDS'
                  Write (6,*) ' User diagnostic duplicated'
                  Write (6,*) ' Anc Ref No ',Anc_Ref_No, &
                              ' Model ',Im_ident
                  ErrorStatus=100
                  CMESSAGE='User diagnostic duplicated'
                  Call Ereport( RoutineName, ErrorStatus, Cmessage )

                End If

              End If

              !     Determine appropriate row number
              If (LModel    == Im_ident   .AND. &
                  LAncRefNo == Anc_Ref_No .AND. &
                  UAncRow   == 0)             THEN

                UAncRow=I    ! Row number found

                !  Over-write pre-existing record
                ANC_record( UAncRow ) = ANC_tmp

                !  Set flag to indicate User record
                ANC_record( UAncRow ) % anc_flag = 'U'

              Else If ( ( LModel > Im_ident  .AND. UAncRow == 0) .OR. &
                        ( LModel == Im_ident .AND.                    &
                          LAncRefNo > Anc_Ref_No  .AND.               &
                          UAncRow == 0)) Then

                UAncRow=I    ! Row number found

                ! This record will be inserted between two pre-existing
                ! records. Create a spare row by moving all subsequent
                ! records up by one row.

                write (6,*) ' User ANCILmaster record to be inserted'
                write (6,*) ' between two pre-existing records.'
                write (6,*) ' Anc Ref No ',Anc_Ref_No,' Model ',Im_ident

                write (6,*) ' Moving records ',UAncRow,':',&
                            ancRecs,' to ',UAncRow+1,':',ancRecs+1

                ! Move records down to create spare row
                Do ID = ancRecs+1, UAncRow+1, -1
                  ANC_record( id ) = ANC_record( id - 1 )
                End Do

                !  Insert new record
                ANC_record( UAncRow ) = ANC_tmp

                !  Set flag to indicate User record
                ANC_record( UAncRow ) % anc_flag = 'U'

              End If

            Else If (LAncRefNo == 0 .AND. UAncRow == 0) Then

              ! This record will be added after all pre-existing records

              write (6,*) ' User ANCILmaster record to be '
              write (6,*) ' inserted after all pre-existing records.'
              Write (6,*) ' Anc Ref No ',Anc_Ref_No,' Model ',Im_ident

              UAncRow = I

              !  Add new record
              ANC_record( UAncRow ) = ANC_tmp

              !  Set flag to indicate User record
              ANC_record( UAncRow ) % anc_flag = 'U'

            End If

          End Do   !  Loop over i

          !  Increment ancRecs if UserANCIL record has been added.
          If (.NOT. OVERWRITE ) Then

            ancRecs = ancRecs + 1

!         Remove eventually
          If (PrintStatus >= PrStatus_Diag .and. mype == 0) then
            write (6,*) ' ancRecs increased to ',ancRecs
            write (6,*) ' Ancil Ref Nos : ', &
                        ANC_record(:) % ancil_ref_number
          endif

          End If  ! If (.not.overwrite)

      End If  !  if .not. UsrAnc

    End If      ! im_ident /=1

  End If      ! Char == '1'

End Do !  Loop over records until im_ident = -1 reached


! Close relevant file and get rid of unit
If ( mype == 0 ) Then
   Close( nftancilm )
   Call Rcf_Free_Unit( nftancilm )
End If

Return
End Subroutine GetAnc_Fields

End Module GetAnc_Fields_Mod
