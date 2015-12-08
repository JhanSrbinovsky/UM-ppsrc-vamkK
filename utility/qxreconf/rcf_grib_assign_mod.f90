! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Assign 'GRIB record' to the correct linked list

Module Rcf_Grib_Assign_Mod

! SUBROUTINE Rcf_Grib_Assign  - File 'current' data into relevant list
!
! Description: This routine files the information in derived type
!              'Current' into one of the predefined 'lists' used for
!              each parameter. (as specified in rcf_grib_lookups.F90)
!
! Method:  Using the center ID no. (Block_1, octet 1) compare the
!          parameter ID of the data within the lookup table.
!          Use the lookup table to find the stash code with which to
!          file the data in the correct 'list'
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/06/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
SUBROUTINE Rcf_Grib_Assign(Current, Lists)

Use Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params

Use Rcf_GRIB_Block_Params_Mod          ! provides type def for Current

Use EReport_Mod, Only :     &
    EReport

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

USE UM_ParVars, Only : &
    mype

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer       :: Current       !\ Pointer to current
                                                  !/ grib record
Type (List_Marker),Intent(InOut) :: Lists(0:grib_max_fields)

! Local variables

Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Assign'
Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport
Integer                          :: I             ! Loop counter
Integer                          :: LookupCol     ! Column in 'Table A'
                                                  ! of rcf_grib_lookups
Integer                          :: ListNo        ! List to store data

Logical                          :: L_Error       ! Error Flag

!=======================================================================
!  Main Routine
!=======================================================================

!=======================================================================
!  Test record in order to assign to correct list
!=======================================================================

! First - Check which originating center data came from
!         And set the lookup column accordingingly
Select Case (Current % Block_1(p_Orig_cntr))

  Case (GrbOrigECMWF)   ! Grib came from ECMWF
    LookupCol = p_ECMWF_IDCol

  Case (GrbOrigUKMO)   ! Grib came from UKMO
    LookupCol = p_UKMO_IDCol

  Case Default          ! Unidentified center code
    Write (Cmessage(1),'(A,I3,1X,A)')                                 &
                             'Originating Center code number :',      &
                             Current % Block_1(p_Orig_cntr)           &
                             ,'not recognised.'
    ErrorStatus = 10
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )

End Select

! Second - Use the lookup to get the Stash Code and List number
ListNo = 0                     ! The default list for 'misc' fields
Current % StashCode = -1
Current % Desc      = "Unknown Parameter   "

If ( PrintStatus >= PrStatus_Diag .And. mype == 0 ) Then
  Write (6,'(A,I3)') "Looking for Parameter :",                       &
               Current % Block_1(p_Param_ID)
End If

L_Error = .True.
File: Do I = 1, p_Max_Rows
  If (Current % Block_1(p_Param_ID) ==                                &
              Param_ID_CrossRef( I) % CrossRefIDs(LookupCol) .AND.    &
      Current % Block_0(p_Tbl_Vers_No) ==                             &
              Param_ID_CrossRef( I) % Table_No(LookupCol) ) Then

    Current % StashCode =                                             &
                  Param_ID_CrossRef(I) % CrossRefIDs(p_STASH_IDCol)
    ListNo = Param_ID_CrossRef( I) % List_No
    Current % Desc = Param_ID_CrossRef(I) % DescText
    L_Error = .False.

    If ( PrintStatus >= PrStatus_Oper .And. mype == 0) Then
      Write (6,'(2A)') "Found Parameter :",                           &
                                     Param_ID_CrossRef(I) % DescText
    End If
    Exit File                       ! found what I want - exit loop

  End If
End Do File

! double check that new record was found in lookup table.
If ( L_Error ) Then
  If ( PrintStatus >= PrStatus_Normal .And. mype == 0 ) Then
    Write (Cmessage(1),'(A,I3,A)') 'Parameter ',                      &
                         Current % Block_1(p_Param_ID),               &
                         ' in GRIB file not found in lookup table'
    ErrorStatus = -10    ! Warning only - unknown data filed in misc
    Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  End If
End If

!=======================================================================
!  Now put field in the selected list
!=======================================================================

Current % Prev   => Lists(ListNo) % End
                               ! Point Prev pointer at end of
                               ! current list.(Null if first entry)
If (Associated(Lists(ListNo) % End)) Then
                               ! If current end of list is a
                               ! valid record (Not first entry)
  Lists(ListNo) % End % Next  => Current
                               ! Point 'next' for previous entry
                               ! at current entry

Else                           ! Else : must be 1st entry
  Lists(ListNo) % Begin  => Current
                               ! Point begining of List at Current
End If

Lists(ListNo) % End      => Current
                               ! Point End of List at (now complete)
                               ! Current Entry
Nullify(Current % Next)        ! Ensure 'Next' is not associated

Lists(ListNo) % LstCount = Lists(ListNo) % LstCount + 1
                               ! Add one to count of list size

Return

End Subroutine Rcf_Grib_Assign
End Module Rcf_Grib_Assign_Mod
