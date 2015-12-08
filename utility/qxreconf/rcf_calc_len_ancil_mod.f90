! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate work space required for ancillary data.

MODULE Rcf_calc_len_ancil_Mod

!  Subroutine Rcf_Calc_len_ancil  - Calculate work space for anc. data.
!
! Description:
!    Determines the work space required for the ancillary data.
!
! Method:
!    For each ancillary fields to be read in (SOURCE=2 in ITEMS
!    namelist), determine length of data to be read in and accumulate
!    to return overall workspace required.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE Rcf_Calc_len_ancil (P_Field, R_Field, P_Rows, Len_Ancil)

USE Submodel_Mod, Only :     &
    N_Internal_Model,        &
    Submodel_Ident

USE Rcf_NRecon_Mod, Only : &
    ReconDatList,              &
    Recondat_Node

USE Rcf_Model_Mod, Only : &
    ZonAvOzone

USE Rcf_Stashcodes_Mod, Only : &
    stashcode_ozone,           &
    stashcode_riv_sequence,    &
    stashcode_riv_direction,   &
    stashcode_riv_storage,     &
    stashcode_prog_sec

USE Ancil_Mod, Only : &
    AncRecs,          &
    Anc_Record

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Diag

USE Ereport_Mod, Only :   &
    Ereport

USE UM_ParVars, Only : &
    mype

IMPLICIT NONE

! Arguments
Integer  :: P_Field         !  Length of field
Integer  :: R_Field         !  Length of River Routing field
Integer  :: P_Rows          !  No of rows
Integer  :: Len_Ancil       !  Total workspace for anc data

! Local variables
Integer  :: i,irec,j1       !  Loop indices
Integer  :: N_Levs          !  No of levels
Integer  :: N_Pseudo_Levs   !  No of pseudo levels
Integer  :: Len_anc_data    !  Length of anc data
Integer  :: StashCode       !  Stash Code
Integer  :: SectionCode     !  Section Code
Integer  :: Sec_Item        !  1000*Section + Stashcode
Integer  :: ErrorStatus     !  Error return code

Character (Len=80) :: CMessage      ! Error return message
Character (Len=*), Parameter :: RoutineName='Rcf_Calc_Len_Ancil'

ErrorStatus = 0
CMessage = ' '

Len_Ancil = 0

I = Submodel_Ident

Do irec = 1, ancRecs

  If (anc_record(irec) % anc_field_read == 1) Then  !  Read Anc field

    len_anc_data = 0
    StashCode    = anc_record(irec) % item_number
    SectionCode  = anc_record(irec) % section_number
    Sec_Item = 1000*SectionCode + StashCode
    Recondat_Node => RecondatList(I,SectionCode)
    ! Recondat is ordered by section and contains a singly link list to each
    ! item number (in order) so only item number less than whats in list needs
    ! to be compared.  This has room for improvement. 
    Do While ( Associated(Recondat_Node % Recondat_Info) .AND.        &
               Recondat_Node % Recondat_Info % Sec_Item .lt. Sec_Item )
      Recondat_Node => Recondat_Node % Next
    End do
    
    If (Recondat_Node % Recondat_Info % Sec_Item == Sec_Item) Then
      N_Levs = Recondat_Node % Recondat_Info % RLevs
      N_Pseudo_Levs = Recondat_Node % Recondat_Info % RPLevs
    Else
      ErrorStatus=1
      Write (CMessage, '(A, I7)') &
      'StashCode is not a valid prognostic variable, Stash Code ',Sec_Item
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    !       Ancillary data is nearly always read in on full grid.
    !       Special cases are as follows :

    Select Case ( SectionCode )

    Case (stashcode_prog_sec)

      Select Case( StashCode )

      Case (stashcode_ozone)
        !If using Zonal Ozone then size is P_Rows not P_Field.
        If ( ZonAvOzone ) Then
          Len_anc_data = N_Levs * N_Pseudo_Levs * P_Rows
        Else
          Len_anc_data = N_Levs * N_Pseudo_Levs * P_Field
        End If

      Case ( stashcode_riv_sequence,      &
        stashcode_riv_direction,     &
        stashcode_riv_storage )
        ! These fields are on an R_Field sized grid.
        Len_anc_data = N_Levs * N_Pseudo_Levs * R_Field
  
      Case Default
        Len_anc_data = N_Levs * N_Pseudo_Levs * P_Field

      End Select

    Case Default
      Len_anc_data = N_Levs * N_Pseudo_Levs * P_Field

    End Select

    If ( Len_anc_data == 0 ) Then  !  Prognostic not in output dump
      ErrorStatus=10
      Write (CMessage, '(A, I7, A)') &
      ' Ancillary prognostic, Stash Code ',Sec_Item, &
      ' not found in output dump.'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    Len_Ancil = Len_Ancil + Len_anc_data

    If (PrintStatus > PrStatus_Diag .and. mype == 0) Then
      write (6,*)
      write (6,*) ' section_code  ',SectionCode
      write (6,*) ' item_code     ',StashCode
      write (6,*) ' n_levs        ',n_levs
      write (6,*) ' n_pseudo_levs ',n_pseudo_levs
      write (6,*) ' len_anc_data  ',len_anc_data
      write (6,*) ' len_ancil     ',len_ancil
    End If

  End If  ! If ancillary field to be read
End Do    ! irec

Return
End Subroutine Rcf_calc_len_ancil

End Module Rcf_calc_len_ancil_Mod
