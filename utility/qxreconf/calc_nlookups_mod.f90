! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate number of lookup entries in ancillary files

Module Calc_nlookups_Mod
IMPLICIT NONE

!  Subroutine Calc_nlookups  - Get no of lookups in ancillary files.
!
! Description:
!    Calculate number of lookup entries (NLOOKUPS) in ancillary files
!
! Method:
!    Loop over ancillary files to be opened and count the number of
!    lookup entries from the fixed headers.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Calc_nlookups (nlookups)

Use Ancil_mod, Only :     &
    anc_record,           &
    ancrecs,              &
    ancfiles,             &
    anc_file,             &
    ancf_unitno

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_UMhead_Mod, Only : &
    LenFixHd

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParVars, Only : &
    mype

Use Rcf_FortranIO_Mod, Only : &
    Max_Filename_Len
USE IO
IMPLICIT NONE

! Arguments
Integer, Intent(Out)  :: nlookups  ! No of lookup entries

! Local variables
Integer       :: i               ! Loop index
Integer       :: j               ! Loop index
Integer       :: icode           ! Return code
Integer       :: ErrorStatus     ! Error return code
Integer       :: IOStatus        ! Return code from 'Inquire'
Integer       :: len_anc_env_var ! Length of env variable

Integer, dimension(:), allocatable :: fixhd ! Fixed Header

Character (Len=80)           :: Cmessage    !
Character (Len=*), Parameter :: RoutineName = 'Calc_nlookups'

Character (Len=Max_Filename_Len) ::  AncFileName  ! Ancillary file name

Logical            :: l_exist       ! T if Ancillary File exists.
Logical            :: l_files_found ! F is an ancillary file missing.


!-----------------------------------------------------------------
! Determine which ancillary files need to be opened
!-----------------------------------------------------------------

anc_file(:) % anc_file_open = 0

do i=1,ancFiles
  do j=1,ancRecs

    if (anc_record(j) % anc_field_read  == 1   .and. &
        anc_record(j) % anc_file_number == anc_file(i) % anc_file_number ) then

       anc_file (i) % anc_file_open = 1
       exit

     endif

  enddo !  j
enddo  !  i

!-----------------------------------------------------------------
! Print out which ancillary files need to be opened
!-----------------------------------------------------------------

If ( PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
  write (6,*) ' '
  write (6,*) 'Ancillary Files to be opened : '
  do i=1,ancFiles  !  Loop over all ancillary files
    if (anc_file(i) % anc_file_open == 1 ) then
      write (6,'(A,I4,1X,A)') &
        'File No ',anc_file(I) % anc_file_number, anc_file (I) % anc_file_title
    endif  !  If file to be opened
  enddo  !  Loop over ancillary files
  write (6,*) ' '
End If

!-----------------------------------------------------------------
! Allocate space for fixed header
!-----------------------------------------------------------------
allocate ( fixhd(LenFixHd) )

!-----------------------------------------------------------------
! Initialise nlookups
!-----------------------------------------------------------------
nlookups = 0
l_files_found = .TRUE.

!-----------------------------------------------------------------
! Open ancillary files, read in fixed header and count number
! of entries in lookup-table
!-----------------------------------------------------------------
do i=1,ancFiles !  Loop over all ancillary files

  if (anc_file(i) % anc_file_open == 1 ) then ! Anc File to be opened

!----------------------------------------
!    Get the file name from env. variable
!----------------------------------------
     Call Fort_Get_Env( anc_file(i) % anc_env_var, 8, AncFileName,  &
                        Max_Filename_Len, icode)

     If ( icode /= 0 ) Then
       !  If env var not found, return code is -1
       icode = abs (icode)
       write (cmessage, '(2A)')                                 &
       'Cannot get ancillary file name from Env Var ',          &
       anc_file(i) % anc_env_var
       Call Ereport ( RoutineName, icode, cmessage)
     Endif

!------------------------------------
!    Check that ancillary file exists
!------------------------------------
     If (mype == 0) Then
       Inquire ( file=AncFileName, exist=l_exist, iostat=IOStatus )

       If ( .not. l_exist .OR. IOStatus /= 0) then
         Write (6,*) 'Ancillary File does not exist.'
         Write (6,*) 'File : ',AncFileName
         Write (6,*) 'File No ',I,' ',anc_file (I) % anc_file_title
         Write (6,*) 'Got ancillary file name from Env Var ',         &
                                            anc_file(i) % anc_env_var
         l_files_found = .FALSE.
         Cycle
       End If
     End If

!----------------------
!    Get len of env var
!----------------------
     len_anc_env_var = len_trim ( anc_file(i) % anc_env_var )

!---------------------------
!    Open the ancillary file
!---------------------------

     Call File_Open (AncF_UnitNo, anc_file(i) % anc_env_var,      &
                     len_anc_env_var, ioNameInEnv, ioOpenReadOnly,&
                     icode)

     if (icode /= 0) then
       Write (6,*) 'Problem opening Ancillary File.'
       Write (6,*) 'File : ',AncFileName
       Write (6,*) 'File No ',I,' ',anc_file (I) % anc_file_title
       Write (6,*) 'Got ancillary file name from Env Var ',         &
                                          anc_file(i) % anc_env_var
       Cmessage = 'Problem opening Ancillary file.'
       Call Ereport ( RoutineName, icode, Cmessage )
     endif

!---------------------------
!    Ensure file is at start
!---------------------------

     call Setpos (AncF_UnitNo, 0, icode)

     if (icode /= 0) then
       Write (6,*) 'Problem in SETPOS for Ancillary File.'
       Write (6,*) 'File : ',AncFileName
       Write (6,*) 'File No ',I,' ',anc_file (I) % anc_file_title
       Write (6,*) 'Got ancillary file name from Env Var ',         &
                                          anc_file(i) % anc_env_var
       Cmessage = 'Error in SETPOS for Ancillary File.'
       Call Ereport ( RoutineName, icode, Cmessage )
     endif

!-----------------------------------
!    Read in the fixed length header
!-----------------------------------
! DEPENDS ON: read_flh
     call Read_flh (AncF_UnitNo, fixhd, LenFixHd, icode, cmessage)

     if (icode /= 0) then
       Write (6,*) 'Problem in reading fixed header from Anc File.'
       Write (6,*) 'File : ',AncFileName
       Write (6,*) 'File No ',I,' ',anc_file (I) % anc_file_title
       Write (6,*) 'Got ancillary file name from Env Var ',         &
                                          anc_file(i) % anc_env_var
       Cmessage='Problem with reading fixed header from Ancillary File.'
       Call Ereport ( RoutineName, icode, Cmessage )
     endif

     Call File_Close (AncF_UnitNo, anc_file(i) % anc_env_var,    &
                     len_anc_env_var, ioNameInEnv)

!------------------------------------------------------
!    Accumulate no of lookup entries in ancillary files
!------------------------------------------------------
     nlookups = nlookups + fixhd (152)

  endif  !  If file to be opened
enddo  !  Loop over ancillary files

!---------------------------------------------------------
! If a file has not been found, generate an error
!---------------------------------------------------------
If ( .NOT. l_files_found ) Then
  ErrorStatus = 10
  Cmessage = 'Ancillary files have not been found - Check output '//&
        'for details'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-----------------------------------------------------------------
! De-allocate space for fixed header
!-----------------------------------------------------------------
deallocate ( fixhd )

Return
End Subroutine calc_nlookups

End Module Calc_nlookups_Mod
