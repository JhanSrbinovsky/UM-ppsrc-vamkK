! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Interface to STASH processing

Module Rcf_Stash_Proc_Mod

!  Subroutine Rcf_Stash_Proc - STASH processing
!
! Description:
! This routine is an interface to most of the stash processing -
! reading in the  stashmaster files and setting up the addressing
! of the output dump
!
! Method:
!   Does some initialisation, then calls rcf_getppx, submodl and
!   rcf_address.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains


Subroutine Rcf_Stash_Proc( EnvDir, UsmEnvDir, stash_in_arg )

Use Ereport_Mod, Only:   &
    Ereport

Use Rcf_Getppx_Mod, Only    :     &
    Rcf_Getppx

Use Submodel_Mod, Only :      &
    a_im,                     &
    n_submodel_partition_max, &
    n_internal_model_max,     &
    Internal_Model_Index

Use Rcf_Model_Mod, Only :         &
    H_Vers

Use Rcf_NRecon_Mod, Only :    &
    DumpProgLevs,             &
    PrimDataLen

Use Rcf_PPX_Info_Mod, Only : &
    STM_Record

IMPLICIT NONE

! Subroutine arguments
Character (Len=*), Intent(In)   :: EnvDir
                                    ! Env. Variable for STASHmaster
                                    ! directory
Character (Len=*), Intent(In)   :: UsmEnvDir
                                    ! Env. Var for user STASHmaster
                                    ! files
Logical, Optional, Intent(In)   :: stash_in_arg
                                    ! T if input STASHmaster

! Local Data
Character (Len=*), Parameter    :: RoutineName = 'Rcf_StashProc'
Character (Len=80)              :: CMESSAGE     ! Error return message
Integer                         :: ErrorStatus  ! +ve = fatal error
Logical                         :: stash_in     ! Local value of stash_in_arg

!--------------------------------------------------------------------

! Set default for STASHMSTR_IN flag
If ( Present (stash_in_arg) ) Then
  stash_in = stash_in_arg
Else
  stash_in = .FALSE.
End If

If ( .NOT. stash_in ) Then

! Initialisation
  ErrorStatus =0

! Initialisation of data length arrays
  H_VERS( :, : )  = 0

  DumpProgLevs(:) = 0
  PrimDataLen (:) = 0

! Nullify STM_Record array
  Nullify( STM_Record )
End If

! Read STASHmaster files into look-up arrays PPXI, PPXC
If (INTERNAL_MODEL_INDEX(A_IM) > 0) Then
  Call Rcf_Getppx( 'STASHmaster_A', EnvDir, stash_in_arg = stash_in )
End If

!Read user STASHmaster files (which may be empty)
Call Rcf_Getppx('             ', UsmEnvDir, .TRUE., stash_in )

Return
End Subroutine Rcf_Stash_Proc

End Module Rcf_Stash_Proc_Mod
