! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Diagnostic tools used with GRIB data

Module Rcf_Grib_Debug_Tools_Mod
!
! Description:
!   Small selection of routines allowing me to print information from
!   the linked lists set up to handle incoming GRIB data
!
! *** NOTE *** - Routines containrd herin have unprotected write
!                statements, therefore calls to this routine should be
!                protected by normal procedures for write statements
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/06/02   Original code. Roddy Sharp (frtz)
!
!-----------------------------------------------------------------------
!  Variables to do specifically with the GRIB record
!-----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

Use Rcf_GRIB_Block_Params_Mod

Use Rcf_GRIB_FldSort_Mod, Only  :  &
  Rcf_GRIB_FldSort

Use Rcf_StashCodes_Mod

Use Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params

Type (Grib_Record),Pointer       :: Current       !\ Pointer to current
                                                  !/ grib record

! Function definitions
!=======================================================================
! *** Actual Subroutines Involved                                  ***
!=======================================================================

contains
!=======================================================================
!***********************************************************************
!* Print some basic information on the entries in                      *
!* each list (except the 'unknown parameters' one                      *
!***********************************************************************
!=======================================================================

SUBROUTINE Grib_Debug_Print_Basics(Lists)

Implicit None

!An array of pointer pairs to form the head and tail of the field lists
Type (List_Marker)               :: Lists(0:grib_max_fields)

Integer                          :: I

!=======================================================================
! Loop through lists displaying basic info              ! debug
!=======================================================================

Do I = 1, grib_max_fields                    ! Loop through all lists

If (.NOT.Associated(Lists(I)%Begin)) Then    ! check list has a start

  Write (6,*) " List is Unassociated and contains no entrys"

Else                                        ! List _has_ a start entry

  Write (6,'(3A,I2,A)') &
             "List ", Trim(Lists(I) % Begin % Desc) , " contains ",  &
              Lists( I) % LstCount, " members"

  Current => Lists(I) % Begin
  Write (6,'(A15,1x,":")',Advance='no') "Parameter ID is"
  Do While (Associated(Current))
    Write (6,'(2x,I3)',Advance='no') Current % Block_1 ( p_Param_ID )
    Current => Current % Next
  End Do   ! while associated current
  Write (6,*)

  Current => Lists(I) % Begin
  Write (6,'(A15,1x,":")',Advance='no') "Type of Levl is"
  Do While (Associated(Current))
    Write (6,'(2x,I3)',Advance='no') Current % Block_1 ( p_Lvl_Type )
    Current => Current % Next
  End Do   ! while associated current
  Write (6,*)

  Current => Lists(I) % Begin
  Write (6,'(A15,1x,":")',Advance='no') "1st Lvl Desc is"
  Do While (Associated(Current))
    Write (6,'(1x,I4)',Advance='no') Current % Block_1(p_Lvl_Desc_1)
    Current => Current % Next
  End Do   ! while associated current
  Write (6,*)

  Current => Lists(I) % Begin
  Write (6,'(A15,1x,":")',Advance='no') "2nd Lvl Desc is"
  Do While (Associated(Current))
    Write (6,'(1x,I4)',Advance='no') Current % Block_1(p_Lvl_Desc_2)
    Current => Current % Next
  End Do   ! while associated current
  Write (6,*)

  Current => Lists(I) % Begin
  Write (6,'(A15,1x,":")',Advance='no') "Start point is "
  Do While (Associated(Current))
    Write (6,'(1x,I7)',Advance='no') Current % Start_pos
    Current => Current % Next
  End Do   ! while associated current
  Write (6,*)

  Current => Lists(I) % Begin
  Write (6,'(A15,1x,":")',Advance='no') "Record Length  "
  Do While (Associated(Current))
    Write (6,'(1x,I6)',Advance='no') Current % Block_0(p_Mes_Len)
    Current => Current % Next
  End Do   ! while associated current
  Write (6,*)

  Write (6,*)

End If ! associated list begin
End Do ! loop over all lists

Return

End SUBROUTINE Grib_Debug_Print_Basics

!=======================================================================
!***********************************************************************
! Print a count of entries for each list                ! debug
!***********************************************************************
!=======================================================================

SUBROUTINE Grib_Debug_ListCounts(Lists)

Implicit None

!An array of pointer pairs to form the head and tail of the field lists
Type (List_Marker)               :: Lists(0:grib_max_fields)

Integer                          :: I

Do I = 1, grib_max_fields
  If (Associated( Lists( I) % Begin )) Then
    Write (6,'(A20,A,I3,A)')                                      &
              AdjustR(Lists( I) % Begin % Desc) , " contains ",   &
              Lists( I) % LstCount, " members"
  Else
    Write (6,'(A20,A,I3,A)')                                      &
              "** No Description **" , " contains ",              &
              Lists( I) % LstCount, " members"
  End If
End Do

Write (6,'(A,I3,A)') "     grib_Misc_field contains ",            &
            Lists(grib_misc_field)% LstCount,  " members"

Return

End SUBROUTINE Grib_Debug_ListCounts

!=======================================================================
!***********************************************************************
!* Print contents of the blocks for a Record passsed in      ! debug   *
!***********************************************************************
!=======================================================================

SUBROUTINE Grib_Debug_Print_Blocks(Record)

Implicit None

Type (Grib_Record), Pointer      :: Record

Integer                          :: I
Character(len=50)                :: cFormat1,cFormat2


cFormat1 = "(A,I2,A,I10)"
cFormat2 = "(A)"

! Loop through the blocks of the entry passed in.
If (Associated(Record)) Then
  ! Block_0
  Write (6,cFormat2) "** Block_0 **"
  Do I = 1,Len_Block_0
    Write (6,cFormat1) "Octet ",I," reads ", Record % Block_0(I)
  End Do

  ! Block_1
  Write (6,cFormat2) "** Block_1 **  -- first 50 elements only"
  Do I = 1,50    ! The rest is _normally_ blank
    Write (6,cFormat1) "Octet ",I," reads ", Record % Block_1(I)
  End Do

  ! Block_2
  Write (6,cFormat2) "** Block_2 **"
  Do I = 1,Len_Block_2
    Write (6,cFormat1) "Octet ",I," reads ", Record % Block_2(I)
  End Do

  ! Block_3
  Write (6,cFormat2) "** Block_3 **"
  Do I = 1,Len_Block_3
    Write (6,cFormat1) "Octet ",I," reads ", Record % Block_3(I)
  End Do

  ! Block_4
  Write (6,cFormat2) "** Block_4 **"
  Do I = 1,Len_Block_4
    Write (6,cFormat1) "Octet ",I," reads ", Record % Block_4(I)
  End Do

  ! Block_R
  Write (6,cFormat2) "** Block_R **"
  Do I = 1,Len_Block_R
    Write (6,*) "Octet ",I," reads ", Record % Block_R(I)
  End Do

Else  ! record was null on calling

  Write (6,cFormat2) "rcf_Grib_Debug_Print_Blocks was passed an " // &
              "unassociated record"
End If     ! Associated (Record)

Return

End SUBROUTINE Grib_Debug_Print_Blocks

End Module Rcf_Grib_Debug_Tools_Mod
