! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads stashmaster file into an array of records

Module Rcf_Getppx_Mod

!  Subroutine Rcf_Getppx - read in stashmaster file into array
!
! Description:
!   Reads a series of records from a stashmaster file into the
!   internal array data-structure and updating the ppxptr reference
!   pointer as required
!
! Method:
!   Normal and User stashmaster entries dealt with seperately
!   Normal entries are assumed to be increasing in section/item no
!   Whereas user stashmaster requires more work to ensure that
!   records are inserted/replaced at the correct point.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Counters for array position
Integer, Private, Target, SAVE   :: RowNumber    = 0   
Integer, Private, Target, SAVE   :: RowNumber_in = 0

CONTAINS

SUBROUTINE Rcf_Getppx( StmsrNam, EnvDir, UsrSTM, stash_in_arg )

USE Submodel_Mod, Only :          &
    Internal_Model_Index

USE Rcf_FortranIO_Mod, Only :     &
    Rcf_Get_Unit,                 &
    Rcf_Free_Unit,                &
    Max_Filename_Len

USE Ereport_Mod, Only :   &
    Ereport

USE Rcf_Ppx_Info_Mod, Only :  &
    ppxptr,               &
    ppxptr_in,            &
    n_ustash,             &
    n_ustash_in,          &
    nrecs_ustash,         &
    nrecs_ustash_in,      &
    ustsfils,             &
    ustsfils_in,          &
    ppxrecs,              &
    ppxrecs_in,           &
    STM_record_type,      &
    STM_record,           &
    STM_record_in

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

USE Rcf_ReadSTM_Mod, Only :   &
    Rcf_ReadSTM

USE Rcf_BcastSTM_Mod, Only :  &
    Rcf_BcastSTM

USE UM_ParVars, Only :   &
    mype,                 &
    nproc

USE version_mod, ONLY :  &
    NDiagP

IMPLICIT NONE

! Arguments
Character (Len=13), Intent(In)  :: StmsrNam    ! Names of stash master
                                               ! files
Character (Len=*), Intent(In)   :: EnvDir      ! Env variable for
                                               ! directory name
Logical, Optional, Intent(In)   :: UsrSTM      ! Logical controling
                                               ! User Stashmaster
Logical, Optional, Intent(In)   :: stash_in_arg ! T if input STASHmaster

! Local Variables
Character (Len=*), Parameter :: RoutineName='Getppx'
Integer            :: I, ID         ! Loop counters
Integer            :: ErrorStatus   ! Error return code
Integer            :: hashcount
Integer            :: IFIL,IREC     ! Do.
Integer            :: nftppxref     ! Unit number for STASHmaster
Integer            :: nftstmstu     ! Unit number for UserSTASH
Integer            :: IOSTATUS
Integer            :: Im_index      !
Integer            :: Im_ident      !
Integer            :: Section       !
Integer            :: Item          !
Integer            :: LModel  ,DM
Integer            :: LSection,DS
Integer            :: LItem   ,DI
Integer            :: USTrow
Integer            :: FirstBlank    !Used to append Upsm filename to dir
Integer            :: RI            ! Row index
Integer            :: NU_recs       ! No. of records in a user psm file
Integer            :: Info          ! For GCOM
Integer            :: msg           ! For GCOM tag
Integer            :: EnvLen        ! Length of EnvDir character string
Integer, Pointer   :: ppxrecs_ptr   ! Pointer for ppxrecs/_in
Integer, Pointer   :: rownumber_ptr ! Pointer for rownumber/_in
Integer, Pointer   :: nrecs_ustash_ptr  ! Pointer for nrecs_ustash/_in
Integer, Pointer   :: n_ustash_ptr      ! Pointer for n_ustash/_in
Integer, Pointer   :: ppxptr_ptr(:,:,:) ! Pointer for ppxptr/_in

Character (Len=Max_Filename_Len) :: UpsmFile
                                    ! Full pathname for user psm files
Character (Len=Max_Filename_Len) :: STASH_MSTR
                                    ! Do. STASH master files

Character (Len=80) :: CMESSAGE      ! Error return message
Character (Len=1)  :: CHAR1
Character (Len=80), Pointer :: ustsfils_ptr(:) ! Pointer for ustsfils/_in

LOGICAL            :: OVERWRITE     ! T if a system stash master record
                                    !is being overwritten by a user rec
Logical            :: UUsrSTM       ! Controlling UserSTASH
Logical            :: l_exist       ! T if STASHmaster file exists
Logical            :: stash_in      ! Local value of stash_in_arg

Type (STM_record_type) :: STM_tmp   ! Temporary record
Type (STM_record_type), Pointer :: STM_record_ptr(:) ! Pointer for STM_record/_in

!- End of header -------------------------------------------------------

ErrorStatus = 0
NU_recs     = 0
IOStatus    = 0

! Set default for UserSTASH flag
If ( Present( UsrSTM ) ) Then
  UUsrSTM = UsrSTM
Else
  UUsrSTM = .False.
Endif
! Set default for STASHMSTR_IN flag
If ( Present (stash_in_arg) ) Then
  stash_in = stash_in_arg
Else
  stash_in = .FALSE.
End If

! Allocate space for the ppxc, ppxi and ppxptr arrays where required.
! If one needs allocating, they all do.
If ( stash_in ) Then
  If ( .Not. Associated ( STM_record_in ) ) Then
    Allocate( STM_record_in( ppxRecs_in ) )

! Initialise important data
    STM_record_in(:) % RowIndex = 0
  End If
Else
  If ( .Not. Associated ( STM_record ) ) Then
    Allocate( STM_record( ppxRecs ) )

! Initialise important data
    STM_record(:) % RowIndex = 0
  End If
End If

! Set up variables and pointers dependent on input/output stash
If ( stash_in ) Then
  ppxrecs_ptr      => ppxrecs_in
  rownumber_ptr    => rownumber_in
  ppxptr_ptr       => ppxptr_in
  nrecs_ustash_ptr => nrecs_ustash_in
  n_ustash_ptr     => n_ustash_in
  ustsfils_ptr     => ustsfils_in
  STM_record_ptr   => STM_record_in
Else
  ppxrecs_ptr      => ppxrecs
  rownumber_ptr    => rownumber
  ppxptr_ptr       => ppxptr
  nrecs_ustash_ptr => nrecs_ustash
  n_ustash_ptr     => n_ustash
  ustsfils_ptr     => ustsfils
  STM_record_ptr   => STM_record
End If

! Calculate length of environment variable to pass to C routine.
envlen = LEN_TRIM(envdir)

!----------------------------------------------------------------------
! Check that the no. of requested diagnostics does not exceed max
!----------------------------------------------------------------------
IF ( stash_in .AND. ppxRecs_in > NDiagP ) THEN
  Write(6,*) 'ERROR: no. of diags. requested exceeds max'
  Write(6,*) 'ppxRecs_in=',ppxRecs_in,' NDIAGP=',NDiagP
  Errorstatus=103
  Cmessage = 'GETPPX: ppxRecs_in > NDIAGP '

  Call Ereport( RoutineName, ErrorStatus, Cmessage )

ELSE IF ( (.NOT. stash_in) .AND. ppxRecs > NDiagP) THEN
  Write(6,*) 'ERROR: no. of diags. requested exceeds max'
  Write(6,*) 'ppxRecs=',ppxRecs,' NDIAGP=',NDiagP
  Errorstatus=104
  Cmessage = 'GETPPX: ppxRecs > NDIAGP '

  Call Ereport( RoutineName, ErrorStatus, Cmessage )
END IF
!----------------------------------------------------------------------

! Are we dealing with non-user STASH?
If ( .NOT. UUsrSTM ) Then
  !---------------------------------------------------------------------
  !Read in records from STASHmaster for current internal model
  !--------------------------------------------------------------------
  ! Get a file unit

  If ( mype == 0 ) Then
    Call Rcf_Get_Unit( nftppxref )

!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename

    Call Fort_Get_Env( EnvDir, EnvLen, STASH_MSTR, Max_Filename_Len,    &
                       ErrorStatus )

    If ( ErrorStatus /= 0 ) Then
      Cmessage = 'Cannot get STASH directory from Env. Vars'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    Endif

    FirstBlank = Len_Trim( STASH_MSTR ) + 1
    STASH_MSTR = Trim( STASH_MSTR )

    STASH_MSTR(FirstBlank:FirstBlank)='/'
    STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam

    !   Check that the STASHmaster file exists
    INQUIRE ( file=STASH_MSTR, exist=l_exist, iostat=IOStatus )

    If ( .not. l_exist ) then
      Write (6,*) 'STASHmaster File does not exist.'
      Write (6,*) ' FILE=',STASH_MSTR
      ErrorStatus=100
      CMESSAGE=' GETPPX: STASHmaster file does not exist.'

      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    !   Open the STASHmaster file
    Open( Unit=NFTPPXREF, File=STASH_MSTR, Iostat=IOStatus )

    If(IOStatus.NE.0) Then
      Write (6,*) 'ERROR in routine GETPPX'
      Write (6,*) &
       &   'CANNOT OPEN STASHmaster FILE, IOSTATUS=',IOStatus
      Write (6,*) 'UNIT=',NFTPPXREF,' FILE=',STASH_MSTR
      ErrorStatus=101
      CMESSAGE=' GETPPX: ERROR OPENING STASHmaster'

      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End If   ! pe0

  Im_ident = 0
  Do While (Im_ident /= -1)
    If ( mype == 0 ) Then
      Read( NFTPPXREF, '(A1)' ) CHAR1
    End If

    msg = 9001
    Call gc_cbcast( msg, 1, 0, nproc, info, char1 )

    If (CHAR1.EQ.'1') Then
      !Read block of records
      If ( mype == 0 ) Then
        Backspace NFTPPXREF
        Call Rcf_ReadStm( STM_tmp, NFTPPXREF, ErrorStatus, CMESSAGE)
      End If

      Call Rcf_BcastSTM( STM_tmp, 0 )

      Im_ident = STM_tmp % model
      Section  = STM_tmp % section
      Item     = STM_tmp % item
      If (Im_ident /= -1) Then      ! Not end of file
        Im_index= INTERNAL_MODEL_INDEX(Im_ident)
        RowNumber_ptr = RowNumber_ptr + 1   !   Increment row number
        ! Assign value to PPXPTR element corresponding to this record
        ppxptr_ptr(Im_ident,Section,Item) = RowNumber_ptr
        !   Transfer data from ppx record to look-up arrays
        STM_record_ptr( RowNumber_ptr ) = STM_tmp

        ! Set row index - indicates values of model,sec,item for
        ! this row
        STM_record_ptr( RowNumber_ptr ) % RowIndex  =  Im_ident*100000 &
                                                       + Section *1000 &
                                                       + Item
        !   Set flag to indicate record originated from ppxref file
        STM_record_ptr( RowNumber_ptr ) % OriginFlag = 'P'

        If (RowNumber_ptr > ppxRecs_ptr) Then
          ErrorStatus = 10
          Write ( Cmessage , '(A, I8)' ) &
            ' PPXI row number exceeds total no. of ppx records ', &
            RowNumber_ptr

          Call Ereport( RoutineName, ErrorStatus, Cmessage )

        End If
      Else      ! Im_ident /= -1
        If ( mype == 0 ) Then
          Close( nftppxref )
          Call Rcf_Free_Unit( nftppxref )
        End If  ! mype == 0
      End If    ! Im_ident /= -1
    End If      ! Char == '1'
  End Do        ! Im_ident /= -1
Else            ! Must be user Stash...
! ----------------------------------------------------------
! Insert user-defined diagnostics into ppxref look-up arrays
! ----------------------------------------------------------

If (nrecs_ustash_ptr > 0) Then
  ! There are user diagnostic records
  ErrorStatus=0
  IOStatus   =0

  If ( mype == 0 ) Then
    ! Get directory name for Upsm files
    Call Fort_Get_Env( EnvDir, EnvLen, UpsmFile, Max_Filename_Len,      &
                       ErrorStatus )

    If (ErrorStatus /= 0) Then
      Cmessage = 'Cannot get User STASH directory from Env. Var.'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    Endif

    FirstBlank = Len_Trim( UpsmFile ) + 1
    UpsmFile   = Trim( UpsmFile )
  End if   ! mype == 0

    ! Loop over user pre-stash master files
  DO IFIL = 1, n_ustash_ptr
    If ( mype == 0 ) Then
      UpsmFile(FirstBlank  :FirstBlank  )='.'
      UpsmFile(FirstBlank+1:FirstBlank+8)=ustsfils_ptr(ifil)

      ! Get a file unit
      Call Rcf_Get_Unit( nftstmstu )

      !   Check that the file exists.
      INQUIRE ( file=UpsmFile, exist=l_exist, iostat=IOStatus )

      If ( .not. l_exist ) then
        Write (6,*) 'USER STASHmaster File does not exist.'
        Write (6,*) 'UNIT=',NFTSTMSTU,' FILE=',UpsmFile
        ErrorStatus=100
        CMESSAGE=' GETPPX: User STASHmaster file does not exist.'

        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

      !   Open user stash master file
      Open( Unit = NFTSTMSTU, File = UpsmFile, IOSTAT = IOStatus)

      If ( IOStatus /= 0 ) Then
        Write (6,*) 'CANNOT OPEN USER PPXREF FILE.IOSTATUS=', &
                                                     IOStatus
        Write (6,*) 'UNIT=',NFTSTMSTU,' FILE=',UpsmFile
        ErrorStatus=101
        CMESSAGE=' GETPPX: ERROR OPENING USER PPXREF'

        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

      !   Read number of records in this file
      Read( NFTSTMSTU, '(I3)' ) NU_recs
    End If  ! mype == 0

    msg = 9002
    Call gc_ibcast( msg, 1, 0, nproc, info, NU_recs )

    !   Read in records from user pre-stash master file
    Do IREC = 1, NU_recs
      !   Initialise OVERWRITE  switch
      OVERWRITE  = .FALSE.
      hashcount = -1
      char1     = '0'

      If ( mype == 0 ) Then
        Do While ( char1 /= '1' )
          Read( NFTSTMSTU, '(A1)' ) CHAR1
          hashcount=hashcount+1
          If (hashcount > 20) Then
            Errorstatus=100
            Cmessage = 'INCORRECT FORMAT IN USER STASHmaster FILE'
            Write (6,*) 'INCORRECT FORMAT IN USER STASHmaster FILE'
            Write (6,*) 'GAP BETWEEN RECORDS TOO LARGE?'
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
        End Do

        !Read block of records
        Backspace NFTSTMSTU
        Call Rcf_READSTM( STM_tmp, NFTSTMSTU, ErrorStatus,CMESSAGE)
      End if

      Call Rcf_BcastSTM( STM_tmp, 0 )

      Im_ident = STM_tmp % model
      Section  = STM_tmp % section
      Item     = STM_tmp % item

      !   Transfer data from ppx record to look-up arrays
      !   No. of records extracted from STASHmaster file(s)= RowNumber.
      USTrow    =   0

!     Bug corrected by DR 16/6 - been there since 4.1 !
!     Only affects User Records that are inserted after all
!     existing records.

!     Do I = 1, RowNumber  - Incorrect.
      Do I = 1, ppxRecs_ptr
        RI      =   STM_record_ptr( I ) % RowIndex
        !     Determine values of model,section,item for this row
        If (RI > 0 .AND. USTrow == 0) THEN
          LModel  =     RI/100000
          LSection=(RI-(RI/100000)*100000)/1000
          LItem   =(RI-(RI/1000  )*1000  )
          !     Check whether previous item is being overwritten
          If (Im_ident == LModel   .AND. &
              Section  == LSection .AND. &
              Item     == LItem        ) THEN

            If      (STM_record_ptr(I) % OriginFlag == 'P') Then
              OVERWRITE =.TRUE.
              If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
              Write (6,*) 'MESSAGE FROM ROUTINE GETPPX:'
              Write (6,*) &
                  'The following PPXREF record has been overwritten by'
              Write (6,*) &
                  'a record read from a user-STASH master file: '
              Write (6,*) 'Internal Model ',Im_ident, &
                  ' Section ',Section,' Item ',Item
              End If

            Else If (STM_record_ptr(I) % OriginFlag == 'U') Then
              Write (6,*) 'ERROR, GETPPX: '
              Write (6,*) 'User diagnostic duplicated'
              Write (6,*) 'Model,Section,Item ', &
                               Im_ident,Section,Item
              ErrorStatus=100
              CMESSAGE='ERROR,GETPPX:user diag duplicated'
              Call Ereport( RoutineName, ErrorStatus, Cmessage )
            End If
          End If

          !     Determine appropriate row number
          If (LModel   == Im_ident .AND. &
              LSection == Section  .AND. &
              LItem    == Item     .AND. &
              USTrow   == 0)             THEN

            USTrow=I    ! Row number found
            !     This record will overwrite a pre-existing record
            !     Insert new record
            STM_tmp % RowIndex = STM_record_ptr( USTrow ) % RowIndex
            STM_record_ptr( USTrow ) = STM_tmp

            ! Set flag to indicate record originated from user psm file
            STM_record_ptr( USTrow ) % OriginFlag = 'U'

          Else If ( ( LModel > Im_ident .AND. USTrow == 0) .OR.  &
                    ( LModel == Im_ident .AND.                   &
                      LSection > Section .AND. USTrow == 0) .OR. &
                    ( LModel == Im_ident .AND.                   &
                      LSection == Section .AND.                  &
                      LItem > Item        .AND. USTrow == 0)) Then
            USTrow=I    ! Row number found
          !This record will be inserted between two pre-existing records
          !Create spare row - move all subsequent records up by one row

            Do ID = RowNumber_ptr+1, USTrow+1, -1
              STM_record_ptr( id ) = STM_record_ptr( id - 1 )
              RI               = STM_record_ptr( id - 1 ) % RowIndex

              !     Determine values of model,section,item for this row
              DM=     RI/100000
              DS=(RI-(RI/100000)*100000)/1000
              DI=(RI-(RI/1000  )*1000  )

              !     Increment PPXPTR for record moved up
              ppxptr_ptr(DM,DS,DI)=ppxptr_ptr(DM,DS,DI)+1
            End Do

            !     Insert new record
            STM_record_ptr( USTrow ) = STM_tmp

            !     Set row index - indicates model,sec,item for this row
            STM_record_ptr( USTrow ) % RowIndex =  Im_ident*100000 &
                                                   + Section *1000 &
                                                   + Item

            !  Set flag to indicate record originated from user psm file
            STM_record_ptr( USTrow ) % OriginFlag = 'U'

            !     Set PPXPTR for the new record
            ppxptr_ptr(Im_ident,Section,Item)=USTrow

          End If

        Else If (RI == 0 .AND. USTrow == 0) Then
          !     This record will be added after all pre-existing records
          USTrow = I

          !     Add new record
          STM_record_ptr( USTrow ) = STM_tmp

          !     Set row index - indicates model,sec,item for this row
          STM_record_ptr( USTrow ) % RowIndex =  Im_ident*100000 &
                                                 + Section *1000 &
                                                 + Item
          !   Set flag to indicate record originated from user psm file
          STM_record_ptr( USTrow ) % OriginFlag = 'U'

          !     Set PPXPTR for the new record
          ppxptr_ptr(Im_ident,Section,Item)=USTrow
        End If
      End Do

      !    Increment RowNumber as UserSTASH record has been added.
      !    don't increment it if a standard record has been overwritten.
      If (.NOT. OVERWRITE ) Then
        RowNumber_ptr = RowNumber_ptr + 1
      End If
    End Do        ! Loop over IREC recs in upsm file

    ! Close relevant file and get rid of unit
    If ( mype == 0 ) Then
      Close( nftstmstu )
      Call Rcf_Free_Unit( nftstmstu )
    End If
  End Do          ! Loop over user psm files
End If          ! nrecs_ustash_ptr > 0

End If  ! Logical Normal or User STASHmaster

NULLIFY( ppxrecs_ptr )
NULLIFY( rownumber_ptr )
NULLIFY( ppxptr_ptr )
NULLIFY( nrecs_ustash_ptr )
NULLIFY( n_ustash_ptr )
NULLIFY( ustsfils_ptr )
NULLIFY( STM_record_ptr )

Return
End Subroutine Rcf_Getppx

End Module Rcf_Getppx_Mod
