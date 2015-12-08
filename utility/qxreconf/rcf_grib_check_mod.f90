! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Perform cconsistancy checks on GRIB data

Module Rcf_Grib_Check_Mod
IMPLICIT NONE

! SUBROUTINE Rcf_Grib_Check  Perform Basic consistancy checks on GRIB
!                            Data being handled.
!-and-
! FUNCTION Find_Match  - recursively traverse a list to find a match
!
! Description: A Routine where checks are performed on the GRIB data
!              read in to ensure it matches expectations.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     12/06/02   Original code. Roddy Sharp (frtz)
!  6.2     20/06/05   Check for consistency of model level
!                     fields added. Paul Earnshaw (frpe)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Check(Lists)

! uses variables and routines from other modules

Use Rcf_GRIB_Block_Params_Mod, Only :    &
  List_Marker,                           &
  Grib_Record,                           &
  p_Lvl_Type,                            &
  p_Param_ID,                            &
  p_Lvl_Desc_1,                          &
  p_Lvl_Desc_2,                          &
  EID_Log_Surf_Press,                    &
  Tb3_Surface,                           &
  Tb3_Pressure,                          &
  Tb3_Hybrid

Use Rcf_GRIB_Lookups_Mod, Only :    &
  grib_max_fields,             &
  grib_Soil_Temp_field,        &
  grib_Soil_Moist_field

Use EReport_Mod, Only :     &
    EReport

Use Rcf_GRIB_T_n_Pstar_H_Interp_Mod, Only :  &
    ak, bk

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(inOut):>
Type (List_Marker), Intent(InOut) :: Lists(0:grib_max_fields)
!An array of pointer pairs to form the head and tail of the field lists

! Local variables

Type (Grib_Record),Pointer       :: Current, Compare

Character (Len=*), Parameter     :: RoutineName='Rcf_Grib_Check'
Character (Len=80)               :: CMessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport
Integer                          :: I,Count,Ref_List,pltype
Integer                          :: nvc,nlp

Logical                          :: Match, Match_Vert
Logical                          :: L_First_Press, L_First_Hybrid

Real                             :: soildepth

!=======================================================================
! Main Routine.
!=======================================================================

!=======================================================================
! Perform Basic Data Integrity Checks
!=======================================================================
                   ! flag to identify first pressure level encountered
L_First_Press    = .True.
                   ! flag to identify first hybrid level encountered
L_First_Hybrid   = .True.

!=======================================================================
!  Loop across all lists, checking they have valid contents
!=======================================================================
Do I = 1, grib_max_fields

  If (Associated(Lists(I) % Begin) ) Then

!=====================================================
! Special check to see if log(pstar) instead of pstar
!=====================================================
! If log(pstar) is field in surface pressure list then
!  alter level information to make it a surface field,
!  not a model level 1 field. Use special routines
!  later to convert from log(pstar) to pstar.
    If (Lists(I)%Begin%Block_1(p_Param_ID)==EID_Log_Surf_Press) Then
      Write(6,*) "Found log(pstar), modifying level information"
      Lists(I) % Begin % Block_1(p_Lvl_Type) = Tb3_Surface
      Lists(I) % Begin % Block_1(p_Lvl_Desc_1) = 0
      Lists(I) % Begin % Block_1(p_Lvl_Desc_2) = 0
    End If

!=====================================
! Special check for soil level fields
!=====================================
! If soil field then store integer index in p_Lvl_Desc_1
!  and store soil depth in p_Lvl_Desc_2
! Note: Fields already sorted
    If ( I == grib_Soil_Temp_field .or. &
         I == grib_Soil_Moist_field ) Then
      Write(6,*) "Modifying soil level information"
      Current => Lists( I ) % Begin
      count   =  0
      Do While ( Associated( Current ) )
        count = count + 1
        soildepth = Current % Block_1 ( p_Lvl_Desc_2 ) -            &
                    Current % Block_1 ( p_Lvl_Desc_1 )
        Current % Block_1 ( p_Lvl_Desc_1 ) = count
        Current % Block_1 ( p_Lvl_Desc_2 ) = soildepth
        Current => Current % Next
      End Do
    End If

!==========================================
! Check No. 1 - Consistant Vertical Levels
!==========================================
! Check 1.1 Levels not repeated for a parameter
!       1.2 List count matches no. of levels found
!       1.3 List counts for different multi-level params match
!       1.4 Levels for different multi-level params also match

    ! Check to see if list lies on Pressure levels or Model levels
    pltype=Lists(I) % Begin % Block_1(p_Lvl_Type)
    If ( pltype == Tb3_Pressure .or. pltype == Tb3_Hybrid ) Then

      ! Check 1.1 - Levels not repeated for a parameter
      !-------------------
      Current => Lists( I ) % Begin

      Do While ( Associated( Current % Next ))

        Match = Find_Match (Current, Current % Next, p_Lvl_Desc_1)

        If (Match) Then     ! This level has same descriptor as another
          Write (CMessage(1),'(A)') "List " //                        &
              Lists( I ) % Begin % Desc     //                        &
              " :Two levels have the same value"
          ErrorStatus = 10
          Call EReport( RoutineName, ErrorStatus, Cmessage(1))
        End If

        Current => Current % Next         ! move to next item in list I

      End Do

      ! Check 1.2  - count
      !-------------------
      Current => Lists( I ) % Begin
      count   =  0

      Do While ( Associated( Current ))
        count = count + 1
        Current => Current % Next
      End Do

      If (Count /= Lists( I ) % LstCount) Then
          Write (CMessage(1),'(A)') "List " //                        &
              Lists( I ) % Begin % Desc     //                        &
              " :  Count of entries does not match recorded count"
        ErrorStatus = -20
        Call EReport( RoutineName, ErrorStatus, Cmessage(1))
        Lists( I ) % LstCount = Count
      End If

      !Is it the first pressure or model level based list found ?
      If ( ( pltype == Tb3_Pressure .and. L_First_Press) .or. &
           ( pltype == Tb3_Hybrid   .and. L_First_Hybrid) ) Then

        ! Are pressure levels and hybrid levels mixed up in same file?
        If (.not. ( L_First_Hybrid .and. L_First_Press)) Then
          Write (CMessage(1),'(A)') "List " //                       &
              Lists( I ) % Begin % Desc     //                       &
              " :  Cannot mix up pressure and model level fields"
          ErrorStatus = 20
          Call EReport( RoutineName, ErrorStatus, Cmessage(1))
        End If

        ! Check 1.2.1 - Vertical Coordinates consistent across records
        !-------------------
        If ( pltype == Tb3_Hybrid ) Then
          Current => Lists( I ) % Begin
          If ( Associated( Current%Next )) Then
            Match_Vert = Find_Match_Vert (Current % Next, Current)
                  If (.Not. Match_Vert) Then
                  Write (CMessage(1),'(A)') "List " //               &
                      Lists( I ) % Begin % Desc     //               &
                      " : Vertical Coordinates not consistent"
              ErrorStatus = 50
              Call EReport( RoutineName, ErrorStatus, Cmessage(1))
            End If
          End If
        End If

        ! record list as 'reference list'
        Ref_list = I
        If ( pltype == Tb3_Pressure) L_First_Press    = .False.
        If ( pltype == Tb3_Hybrid  ) L_First_Hybrid   = .False.

        ! Set ECMWF level definitions using Vertical Coordinates
        If ( pltype == Tb3_Hybrid ) Then
          Current => Lists( I ) % Begin
          nvc=Current%Num_Vert    !\Number of vertical coordinates
          nlp=nvc/2               !/ nvc = 2*(nlevs+1)
          Allocate(ak(nlp-1),bk(nlp-1))
          ak = 0.5 * ( Current%VertCoords(1:nlp-1) + &
                       Current%VertCoords(2:nlp) )
          bk = 0.5 * ( Current%VertCoords(nlp+1:nvc-1) + &
                       Current%VertCoords(nlp+2:nvc) )
        End If

      Else                          ! make comparisons against ref list

        ! Check 1.3 - Match reference levels
        !-----------------------------------
        Current => Lists( I ) % Begin
        Compare => Lists( Ref_list ) % Begin

        Do While ( Associated( Current ))

          Match = Find_Match (Current, Compare, p_Lvl_Desc_1)

          If (.Not.Match) Then   ! This level did not match one in Ref
          Write (CMessage(1),'(A)') "List " //                        &
              Lists( I ) % Begin % Desc     //                        &
              " :  level did not match one in ref field"
            ErrorStatus = 30
            Call EReport( RoutineName, ErrorStatus, Cmessage(1))
          End If

          Current => Current % Next       ! move to next item in list I
        End Do

        ! Check 1.4 - Compare counts
        !-----------------------------------
        If (Lists( I ) % LstCount /= Lists( Ref_list ) % LstCount) Then
          Write (CMessage(1),'(A)') "List " //                        &
              Lists( I ) % Begin % Desc     //                        &
              " :  level count did not match reference"
          ErrorStatus = 40
          Call EReport( RoutineName, ErrorStatus, Cmessage(1))
        End If

        ! Check 1.4.1 - Vertical Coordinates consistent across records
        !-------------------
        If ( pltype == Tb3_Hybrid ) Then
          Current => Lists( I ) % Begin
          Compare => Lists( Ref_list ) % Begin
          If ( Associated( Current )) Then
            Match_Vert = Find_Match_Vert (Current, Compare)
                  If (.Not. Match_Vert) Then
                  Write (CMessage(1),'(A)') "List " //               &
                      Lists( I ) % Begin % Desc     //               &
                      " : Vertical Coordinates did not match reference"
              ErrorStatus = 60
              Call EReport( RoutineName, ErrorStatus, Cmessage(1))
            End If
          End If
        End If

      End If  ! Check for First level

    End If ! level check

  End If ! Begin of list was Associated
End Do  ! loop over all lists

Return

End Subroutine Rcf_Grib_Check

!-----------------------------------------------------------------------
!  Function Find_Match
!-----------------------------------------------------------------------

Recursive Function Find_Match ( Current, Compare, Criteria )          &
                                   Result ( Match )

! Description: Recursively traverse a list, exiting on positive match
!
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
Use Rcf_GRIB_Block_Params_Mod, Only :    &
  Grib_Record

  USE ereport_mod, ONLY : ereport
IMPLICIT NONE
! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)                      :: Criteria

!< Array  arguments with intent(in):>
Type (Grib_Record), Pointer              :: Current, Compare

! Local variables

Logical                          :: Match

!-----------------------------------------------------------------------
! Begin routine
!-----------------------------------------------------------------------

Match = .False.

If ( Current % Block_1(Criteria) == Compare % Block_1(Criteria) ) Then
  Match = .True.
Else
  If ( Associated ( Compare % Next ) ) Then
    Match = Find_Match ( Current, Compare % Next , Criteria )
  End If
End If


End Function Find_Match

!-----------------------------------------------------------------------
!  Function Find_Match_Vert
!-----------------------------------------------------------------------

Recursive Function Find_Match_Vert ( Current, Compare )          &
                                   Result ( Match )

! Description: Recursively traverse a list, exiting on negative match
!
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  6.2     16/06/05   Original code. Paul Earnshaw (frpe)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
Use Rcf_GRIB_Block_Params_Mod, Only :    &
  Grib_Record

USE ereport_mod, ONLY : ereport
IMPLICIT NONE
! Subroutine arguments

!< Array  arguments with intent(in):>
Type (Grib_Record), Pointer              :: Current, Compare

! Local variables

Logical                          :: Match, temp_Match
Integer                          :: I      ! loop counter

!-----------------------------------------------------------------------
! Begin routine
!-----------------------------------------------------------------------

Match = .True.

If (Current % Num_Vert == Compare % Num_Vert ) Then
  Do I=1,Current % Num_Vert
    If (Current % VertCoords(I) /= Compare % VertCoords(I)) &
           Match = .False.
  End Do
Else
  Match = .False.
End If

If ( Match ) Then
  If ( Associated ( Current % Next ) ) Then
    Match = Find_Match_Vert ( Current % Next, Compare )
  End If
End If

End Function Find_Match_Vert

End Module Rcf_Grib_Check_Mod
