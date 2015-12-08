! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read in and count number of ANCILmaster records.

MODULE hdancilm_mod

!  Subroutine HdAncilM  - Count number of ANCILmaster reords.
!
! Description:
!   Count the number of ANCILfields or ANCILfile records.
!
! Method:
!    For an ANCILmaster file (including User files) read in all
!    records and count number of records. The H4 label in the file
!    indicates whether it is records for ancillary fields or files.
!    The total number of records for fields and files are stored in
!    max_ancRecs and max_ancFiles respectively.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE hdancilm ( ancmsrnam , envdir, usranc )

Use Ancil_Mod, Only :     &
    max_ancRecs,          &
    max_ancFiles

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_FortranIO_Mod, Only :        &
    Rcf_Get_Unit,                    &
    Rcf_Free_Unit,                   &
    Max_Filename_Len

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Min

USE UM_ParVars, Only : &
    mype,               &
    nproc

IMPLICIT NONE

! Arguments
CHARACTER (LEN=13) :: ancmsrnam    ! Name of ANCILmaster file
CHARACTER (LEN=9)  :: envdir       ! Name of Env. Var for directory
                                   ! for ANCILmaster
LOGICAL, OPTIONAL  :: usranc       ! Logical controlling User Ancilmaster

! Local parameters:
INTEGER, PARAMETER :: no_version = -9999 ! Value returned by function
                                         ! get_umversion if environment
                                         ! variable VN not set

! Local Scalars
Integer            :: info
Integer            :: msg
Integer            :: ErrorStatus    ! Error return code
Integer            :: FirstBlank     ! Position of first blank found
Integer            :: IOStatus
Integer            :: nft            ! Fortran unit number
Integer            :: model_no       ! Sub-model number
Integer            :: um_version     ! Version of UM
Integer            :: um_revision    ! Revision of UM
Integer            :: acm_version    ! Version of ANCILmaster file
Integer            :: acm_revision   ! Revision of ANCILmaster file
Integer            :: Len_AncMsrNam  ! Length of ANCILmaster filename
Integer            :: First_value    ! First value read from record

Character (Len=80) :: Cmessage
Character (Len=8)  :: c_um_version   ! UM version as string
Character (Len=1)  :: char1
Character (Len=1)  :: c_model_no     ! ANCILmaster sub-model no
Character (Len=8)  :: c_acm_version  ! ANCILmaster version string

Character (Len=*), Parameter :: RoutineName = 'HdAncilM'

Character (Len=Max_Filename_Len) :: ANCIL_MSTR
                                 ! Full pathname for ANCIL master files

Logical            :: found_version  !Indicates presence of ACM version
Logical            :: l_exist        ! T if ANCILmaster file exists.
Logical            :: L_fields       ! T if ANCILmaster for anc fields
Logical            :: L_files        ! T if ANCILmaster for anc files
Logical            :: found_H1       ! T if H1 header found
Logical            :: found_H4       ! T if H4 header found
Integer            :: N_Recs         ! No of records counted
CHARACTER (LEN=12) :: c_acm_type     ! UM version as string
LOGICAL            :: uusranc        ! Controlling User AncilMaster

! Function & Subroutine calls:
Integer            :: get_um_version

External Fort_Get_Env

! Set default for UserANCIL flag
IF ( PRESENT( usranc ) ) THEN
  uusranc = usranc
ELSE
  uusranc = .FALSE.
ENDIF


!-------------------------------------------------------------------
! Most of the processing here only needs to be done by pe0 as it
! involves file manipulations; the results will be broadcast at the end.
!-------------------------------------------------------------------
pe0: If (mype == 0 ) Then

  IOStatus = 0

  ! Get a unit number for ANCILmaster file
  Call Rcf_Get_Unit( nft )

  ! Get directory name for ANCILmaster file
  Call Fort_Get_Env( EnvDir, 9, ANCIL_MSTR, Max_Filename_Len,   &
                     ErrorStatus )

  If ( ErrorStatus /= 0 ) Then
    cmessage = 'Cannot get ANCIL directory from Environment Vars'
    Call Ereport ( RoutineName, ErrorStatus, cmessage )
  Endif

  ! Append filename
  FirstBlank = Len_Trim( ANCIL_MSTR ) + 1
  Len_AncMsrNam = Len_Trim( AncMsrNam )

  ancil_mstr = TRIM( ancil_mstr )
  IF (uusranc) THEN
    ancil_mstr( firstblank:firstblank ) = '.'
  ELSE
    ancil_mstr( firstblank:firstblank ) = '/'
  END IF

  ancil_mstr( firstblank+1:firstblank+len_ancmsrnam ) = ancmsrnam

  !  Check that ANCILmaster file exists.
  INQUIRE ( file=ANCIL_MSTR, exist=l_exist, iostat=IOStatus )

  If ( .not. l_exist ) then
    Write (6,*) 'ANCILmaster File does not exist.'
    Write (6,*) 'FILE=',ANCIL_MSTR
    ErrorStatus=10
    CMESSAGE=' HDANCILM: ANCILmaster file does not exist.'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  !  Open the ANCILmaster file.
  Open( Unit=nft, File=ANCIL_MSTR, Iostat=IOStatus )
  If (PrintStatus >= PrStatus_Normal ) Then
    Write(6,*) '!!!! ANCIL_MSTR ',ANCIL_MSTR
  End If

  If (IOStatus /= 0) THEN
    CMESSAGE= &
    'Error opening ANCILmaster file, routine HDANCILM'
    If ( PrintStatus >= PrStatus_Min ) Then
      Write (6,*) &
      'HdAncilM: Fortran Error Response = ',IOStatus, &
      ' Opening ANCILmaster file ',AncMsrNam
    End If
    ErrorStatus = 20
    Call Ereport ( RoutineName, ErrorStatus, cmessage )
  Endif

  !    Get the UM version from the environment variable $VN.
! DEPENDS ON: get_um_version
      um_version = get_um_version()

      IF ( um_version .ne. no_version ) THEN
    !     Now check through the header section of the ANCILmaster
    !     file looking for H1, H3 and H4

    found_version = .false.
    found_H1 = .false.
    found_H4 = .false.

    Read( nft, '(A1)') char1

    Do While( char1 .eq. 'H' .or. char1 .eq. '#')

      If (char1 .eq. 'H') Then
        Backspace nft
        Read (nft, '(1X, A1)') char1

        If (char1 == '1') Then   !  H1 found in header

          !  This line starts with H1 and contains the sub-model
          !  number the ANCILmaster file is for.
          !  H1| SUBMODEL_NUMBER=1

          Backspace nft
          Read (nft,'(20x,a1)') c_model_no
          Read (c_model_no,'(I1)') model_no

          found_H1 = model_no == 1 .or. model_no == 2

          if (.not. found_H1) Then  ! Invalid H1 header found ?
            Write (cmessage, '(3A)') 'Invalid Sub-Model No ',c_model_no, &
            ' in ANCILmaster header H1.'
            ErrorStatus = 30
            Call Ereport ( RoutineName, ErrorStatus, cmessage )
          endif

        End If  !  char1 == '1'

        If (char1 == '3') Then   !  H3 found in header

          !  This line starts with H3 and should indicate
          !  the ANCILmaster version. The line should look like
          !  H3| UM_VERSION=4.3

          found_version = .true.
          Backspace nft
          Read (nft, '(15x,a8)') c_acm_version
          Read(c_acm_version, '(i1,1x,i1)' ) &
               acm_version, acm_revision
          acm_version = acm_version*100 + acm_revision

          !     Now perform the check against the UM version

          If (acm_version .ne. um_version) then

            Write (cmessage, '(A)') &
            'HDANCILM : UM version and ANCILmaster version differ'
            If ( PrintStatus >= PrStatus_Min ) Then
              Write (6,*) 'Version of ANCILmaster file (' &
              ,acm_version, &
              ') does not match UM version (' &
              ,um_version,') in file ',AncMsrNam
            End if
            ErrorStatus = -40
            Call Ereport ( RoutineName, ErrorStatus, cmessage )

          End If  ! version check

        End If  ! char1 == '3'

        If (char1 == '4') Then   !  H4 found in header

          !  This line starts with H4 and should indicate the type
          !  of ANCILmaster file. The line should look like
          !  H4| TYPE=ANCIL_FIELDS  or
          !  H4| TYPE=ANCIL_FILES

          Backspace nft
          Read (nft, '(9x,a12)') c_acm_type
          L_fields = c_acm_type .eq. 'ANCIL_FIELDS'
          L_files  = c_acm_type .eq. 'ANCIL_FILES'
          found_H4 = l_fields .or. l_files

          if (.not. found_H4) then ! Invalid H4 header found ?
            Write (cmessage, '(2A)') &
            'Invalid H4 label in ANCILmaster file : ',c_acm_type
            ErrorStatus = 50
            Call Ereport ( RoutineName, ErrorStatus, cmessage )
          endif

        End If  ! char1 == '4'

      End If                 ! char1 == 'H'

      Read (nft, '(a1)') char1

    End Do

    If (.not. found_version) Then

      If ( PrintStatus >= PrStatus_Min ) Then
        Write (6,*) &
        'HDANCILM : No ANCILmaster version available; Unable to'
        Write (6,*) &
        'check against UM version for file ',AncMsrNam
      End If

      cmessage = 'No ANCILmaster version available'
      ErrorStatus = -60
      Call Ereport( RoutineName, ErrorStatus, cmessage )

    End If

    !     For safety, rewind to the start of the ANCIL file.
    Rewind (nft)

  End If

  ! Check that both H1 and H4 headers have been found.
  If (.not. found_H1 ) Then
    Write (Cmessage, '(A)') 'H1 header not found in AncilMaster file.'
    ErrorStatus = 70
    Call Ereport ( RoutineName, ErrorStatus, cmessage )
  End If
  If (.not. found_H4 ) Then
    Write (Cmessage, '(A)') 'H4 header not found in AncilMaster file.'
    ErrorStatus = 80
    Call Ereport ( RoutineName, ErrorStatus, cmessage )
  End If

  ! Count records - N_Recs is counter
  N_Recs = 0
  First_value = 0
  Do While (First_value /= -1 )
    Read( nft, '(A1)' ) char1
    If ( char1 == '1' ) Then   ! Read first value from record
      Backspace nft
      Read( nft, '(2X,I5)' ) First_value
      If ( First_value /= -1 ) Then
        N_Recs = N_Recs + 1
      End If
    End If
  End Do

  if (L_fields) then
    max_ancRecs = max_ancRecs + N_Recs
  endif
  if (L_files) then
    max_ancFiles = max_ancFiles + N_Recs
  endif

  If (PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'HdAncilM : No of ANCILmaster records in ', &
                 AncMsrNam, N_Recs
  End If

  Close( Unit = nft )
  Call Rcf_Free_Unit( nft )

End If pe0

!-------------------------------------------------------------------
! Now to broadcast the the variables calculated
!-------------------------------------------------------------------
msg = 7501
Call gc_ibcast( msg, 1, 0, nproc, info, max_ancRecs )

msg = 7502
Call gc_ibcast( msg, 1, 0, nproc, info, max_ancFiles )

Return
End Subroutine HdAncilM

End Module HdAncilM_Mod
