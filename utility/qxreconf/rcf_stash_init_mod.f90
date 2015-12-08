! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisation of STASH related variables

Module Rcf_Stash_Init_Mod

!  Subroutine Rcf_Stash_Init - initialisation of STASH variables
!
! Description:
!   Reads headers for stashmasters and calls rcf_stash_proc to
!   read stashmaster and set up output dump addressing.
!
! Method:
!   Calls rcf_hdppxrf and rcf_stash_proc
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

Subroutine Rcf_Stash_Init( )

Use Rcf_FortranIO_Mod, Only : &
    Max_Filename_Len

USE UM_ParVars, Only : &
    mype

Use Rcf_hdppxrf_mod, Only : &
    Rcf_hdppxrf

Use Rcf_Ppx_Info_Mod, Only: &
    ppxrecs,                &
    ppxrecs_in

Use Rcf_Stash_Proc_mod, Only : &
    Rcf_Stash_Proc

Use Submodel_Mod, Only : &
    N_Internal_Model,    &
    Internal_Model_List, &
    Atmos_IM

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Implicit None

! Local variables
Integer                      :: i
Integer                      :: icode
Logical                      :: stash_in  ! T if STASHMSTR_IN is available
! STASHmaster used for output dump
Character (len=9)    :: stashmaster_out = 'STASHMSTR'
! STASHmaster used for input dump
Character (len=12)   :: stashmaster_in  = 'STASHMSTR_IN'
! Temporary file to try reading env to.
Character (len=Max_Filename_len) :: temp_name
! Routine name
Character (len=*), Parameter :: routinename = 'rcf_stash_init'

!-------------------------------------------------------------------
! Find whether STASHMSTR_IN is available which is currently optional
!-------------------------------------------------------------------
Call fort_get_env(stashmaster_in, 12, temp_name, max_filename_len, icode)
! Check return code and set stashmaster_in to stashmaster_out if problem
If (icode /= 0) Then
  stashmaster_in = stashmaster_out
  stash_in = .FALSE.
Else
  stash_in = .TRUE.
End If

!-------------------------------------------------------------------
! Read in STASHmaster headers
!-------------------------------------------------------------------
Do i = 1, N_Internal_Model
  If ( Internal_Model_List(i) == Atmos_IM ) Then
    Call Rcf_hdppxrf( 'STASHmaster_A',stashmaster_out )
    Call Rcf_hdppxrf( 'STASHmaster_A',stashmaster_in, stash_in_arg = .TRUE. )
  End If
End Do

! User STASHmaster header
Call Rcf_hdppxrf( '             ','         ', .TRUE. )
Call Rcf_hdppxrf( '             ','         ', .TRUE., .TRUE. )

If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
  Write (6,*) 'Rcf_Stash_Init : Total No of STASHmaster records ', &
               ppxRecs
  If ( stash_in ) Then
! Non-ENDGame users do not need to see this.
    Write (6,*) 'Rcf_Stash_Init : Total No of input STASHmaster records ', &
                 ppxRecs_in
  End If
End If

!---------------------------------------------------------------------
! Read in the STASHmasters, and do the output dump addressing
!---------------------------------------------------------------------
Call Rcf_Stash_Proc( stashmaster_out, 'USTSHMSTR'   )
! For now use USTSHMSTER as well for userSTASH since we havent yet
! dealt with it in rcf_readnl_ustsnum to handle one.
Call Rcf_Stash_Proc( stashmaster_in , 'USTSHMSTR', .TRUE.)

End Subroutine Rcf_Stash_Init

End Module Rcf_Stash_Init_Mod
