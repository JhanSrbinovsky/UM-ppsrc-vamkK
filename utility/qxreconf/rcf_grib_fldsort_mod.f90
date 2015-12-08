! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ order a linked list order using criteria provided

Module Rcf_Grib_FldSort_mod

! SUBROUTINE rcf_Grib_FldSort  Establish what information is contained
!                               within the GRIB file.
!
! Description: Routine to re-order the contents of a linked list of
!              GRIB field headers in either ascending or descendng
!              order of any one parameter in Block_0
!
! Method: Quite simply, Bubblesort. It compares 2 adjacent entries and
!         if they aren't in the correct order (relative to each other)
!         swap all the pointers so they appear in the opposite order
!         within the list. Repeat this process as you traverse through
!         the list. Then go back to the list start and repeat again
!         until no change occured during a complete pass through the
!         list.
!         Not the most efficient method but simple and suitable
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
SUBROUTINE rcf_Grib_FldSort(Markers,Criteria,Order)

Use Rcf_GRIB_Block_Params_Mod

Use Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params

Use Rcf_Recon_Mod, Only :   &
    GRIB                      ! Logical T if input file is GRIB

IMPLICIT NONE

! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)               :: Criteria ! the field being compared

!< Logical arguments with intent(in):>
Logical, Intent(In), Optional     :: Order  ! The order desired
                                            ! .TRUE. = Ascending

!< Array  arguments with intent(InOut):>
Type (List_Marker), Intent(InOut)       :: Markers

! Local variables

Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_FldSort'
Character (Len=80)               :: Cmessage(3)
Integer                          :: ErrorStatus
Integer                          :: list_count
Type(Grib_Record), Pointer       :: Current,Stored

Logical                          :: No_Change_Pass
Logical                          :: Ascending

!-----------------------------------------------------------------------
!  Initialise bits and bobs.
!-----------------------------------------------------------------------

! Deal with Optional dummy argument 'Order'
If (Present(Order)) Then
  Ascending = Order
Else
  Ascending = .True.     ! Default to Ascending
End If

No_Change_Pass = .False.
list_count = 0

!-----------------------------------------------------------------------
!  Do until a full pass through the list results in no change.
!-----------------------------------------------------------------------
Do While (.Not.No_Change_Pass)

  list_count = list_count + 1

  Current => Markers % Begin             ! Move to start of list
  No_Change_Pass = .True.                ! re-set fo 'no change'

!-----------------------------------------------------------------------
!  Loop through the list members
!-----------------------------------------------------------------------
  Do While (Associated(Current%next))    ! Do while there is an entry
                                         ! Which follows current

!-----------------------------------------------------------------------
!  Check the ordering of current and the entry which follows it
!-----------------------------------------------------------------------

    If ((Ascending.And.                                               &
       (Current%Next%Block_1(Criteria) < Current%Block_1(Criteria)))  &
         .OR.                                                         &
       (.NOT.Ascending.And.                                           &
       (Current%Next%Block_1(Criteria) > Current%Block_1(Criteria)))) &
         Then

!-----------------------------------------------------------------------
!  If order is incorect then swap the two entrys
!-----------------------------------------------------------------------

      No_Change_Pass = .False.

!-----------------------------------------------------------------------
!  Remove entry _after_ current and update links from current to
!  the entry which follows that.(if it exists)
!-----------------------------------------------------------------------

      Stored       => Current%Next      ! use Stored as a pointer
      Current%next => Stored%Next       ! Might be Null(last entry)

      If (Associated(Stored%Next)) Then !  more entries
        Stored%Next%Prev  => Current    ! Point back from next to cur
      Else                              ! Stored was last entry
        Markers % End     => Current    ! Current is new last entry
      End If

!-----------------------------------------------------------------------
!  re-insert Stored between Current and Current%last.
!-----------------------------------------------------------------------

      If (Associated(Current%Prev)) Then ! There are entries B4 this 1
        Current%Prev%Next => Stored
      Else                              ! going in as first in list
        Markers % Begin   => Stored
      End If

      Stored%Next       => Current
      Stored%Prev       => Current%Prev ! This could still be Null
      Current%Prev      => Stored

    Else                                ! Pair were correctly ordered
      Current  => Current%Next          ! Move to next entry
    End If

  End Do

End Do


!-----------------------------------------------------------------------
!  Done : - return
!-----------------------------------------------------------------------

Return

End Subroutine rcf_Grib_FldSort
End Module rcf_Grib_FldSort_mod
