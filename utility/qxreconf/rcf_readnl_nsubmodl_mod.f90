! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads the NSUBMODL namelist

Module Rcf_Readnl_Nsubmodl_Mod

!  Subroutine Rcf_Readnl_Nsubmodl - reads the NSUBMODL namelist.
!
! Description:
!  Rcf_Readnl_Nsubmodl initialises the model with information specIfying
!  internal model and submodel partitions for the run, which is
!  required for control of coupling when more than one internal model
!  is present.
!
! Method:
!   Data read into the Submodel_Mod module.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Readnl_Nsubmodl( nft )

Use Submodel_Mod   ! Almost all of it...

Use Ereport_Mod, Only : &
    Ereport

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParVars, Only : &
    mype

IMPLICIT NONE

Integer      ::  Nft     ! unit number for namelist

! Local constants
Character (Len=*), Parameter :: RoutineName = 'Um_Submodel_Init'

! Local scalars:
Integer            ::  &
     s,                & ! submodel loop
     i,                & ! internal model loop
     sm,               & ! submodel identifier
     im,               & ! internal model identifier
     sm_prev,          & ! previous submodel identifier
     im_prev             ! previous internal model identifier

Integer            :: ErrorStatus
Character (Len=80) :: Cmessage

! 1. Initialise lists before obtaining values for this experiment.

Do i = 1, N_Internal_Model_Max
   Internal_Model_List(i)      = 0
   Submodel_For_IM(i)          = 0
End Do   ! i over internal model list

Do im = 1, Internal_ID_Max
   Submodel_Partition_Index(im) = 0
   Internal_Model_Index(im)     = 0
   Last_IM_in_SM(im)            = .False.
EndDo   ! im over internal model ids

Do s = 1, N_Submodel_Partition_Max
   Submodel_Partition_List(s) = 0
   Submodel_For_SM(s)         = 0
EndDo  ! s over submodel list

Do sm = 1, Submodel_ID_Max
   N_Internal_For_SM(sm)      = 0
EndDo  ! sm over submodel ids

! 2. Obtain internal model and submodel identIfiers from umui
!    generated namelist.

Read( nft, NSUBMODL)
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  Write(6, Nml=nsubmodl )
End If

! 3. Check umui supplied values
!
! 3.1 Check for umui supplied dimensions against parameter maxima.

If ( N_Internal_Model > N_Internal_Model_Max ) Then
   Write (Cmessage, '(A, I3, A)' )                                   &
        'Too many internal models: ', N_Internal_Model,              &
        '. You need to increase N_Internal_Model_Max'
   ErrorStatus=1       ! Set error flag

   Call Ereport( RoutineName, ErrorStatus, Cmessage )
EndIf

! 3.2 Check umui suppiled values are valid

Do i = 1, N_Internal_Model ! loop over internal models

  im = Internal_Model_List(i) ! internal model identIfier
  If (im <= 0 .OR. im > Internal_ID_Max ) Then
   Write(Cmessage, '(A, I3, A)')                                     &
        'Illegal internal model identifier: ',im,                    &
        '. Check values in namelist NSUBMODL'
   ErrorStatus=1       ! Set error flag

   Call Ereport( RoutineName, ErrorStatus, Cmessage )
  EndIf

  sm = Submodel_For_IM(i)     ! submodel for this internal model
  If (sm <= 0 .OR. sm > Submodel_ID_Max ) Then
   Write(Cmessage, '(A, I3, A)')                                     &
        'Illegal submodel dump identifier: ',sm,                     &
        '. Check values in namelist NSUBMODL'
   ErrorStatus=1       ! Set error flag

   Call Ereport( RoutineName, ErrorStatus, Cmessage )
  EndIf

EndDo ! i=1,N_Internal_Model

! 4. Form internal model and submodel description arrays.

sm_prev = 0             ! Null value of submodel identIfier
N_Submodel_Partition=0  ! Count no. of submodel partitions

Do i = 1, N_Internal_Model ! loop over internal models

  im = Internal_Model_List(i)   ! internal model identIfier
  sm = Submodel_For_IM(i)       ! submodel for this internal model
  Internal_Model_Index(im) = i  ! sequence no. for STASH arrays

  If (sm /= sm_prev) Then  ! new submodel

     N_Submodel_Partition = N_Submodel_Partition+1
     Submodel_Partition_List(N_Submodel_Partition) = sm

    ! Since this is a new submodel, the previous internal model must be
    ! the last internal model in its submodel partition.
     If (N_Submodel_Partition > 1) Then ! Not first dump
        Last_IM_in_SM(im_prev) = .True.
     End If

  EndIf                   ! test on new submodel
  Submodel_For_SM(IM) = N_Submodel_Partition

  Submodel_Partition_Index(im) = sm
  N_Internal_For_SM(sm) = N_Internal_For_SM(sm) + 1

  im_prev = im
  sm_prev = sm

EndDo ! i=1,N_Internal_Model

Last_IM_in_SM(im) = .True.  ! last im in list is last im in sm

! 5. Check calculated dimensions against parameter maxima.

If (N_Submodel_Partition > N_Submodel_Partition_Max ) Then
   Write(Cmessage, '(A, I3, A)')                                     &
     'Too many submodels: ', N_Submodel_Partition,                   &
     '. You need to increase N_Submodel_Partition_Max'
   ErrorStatus=1       ! Set error flag

   Call Ereport( RoutineName, ErrorStatus, Cmessage )
EndIf

Return

End Subroutine Rcf_Readnl_Nsubmodl
End Module Rcf_Readnl_Nsubmodl_Mod
