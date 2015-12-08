! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read in the input grid land-sea masks

Module Rcf_ReadLSMIn_Mod
IMPLICIT NONE

!  Subroutine Rcf_ReadLSMin - read land-sea masks
!
! Description:
! This module is designed to read-in and setup an input land-sea
! mask for multiple PEs. IE, both local and global LSMs are read
! and appropriate sizes calculated. It is assumed that space has
! already been allocated for them, but this is checked.
!
! Method:
!  Local LSM is read in on each PE. This is gathered on PE 0 and
!  broadcast so that all pes have a copy of the global LSM.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!

Contains

Subroutine Rcf_ReadLSMIn(Dump_hdr, fields, field_count )

Use Rcf_Lsm_Mod, Only : &
    glob_lsm_in,    &
    local_lsm_in,   &
    local_land_in,  &
    glob_land_in

Use Rcf_UMhead_Mod, Only : &
    UM_header_type

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars, Only : &
    mype,               &
    nproc,              &
    lasize,             &
    glsizep,            &
    gc_all_proc_group,  &
    current_decomp_type

USE decomp_params, ONLY : &
    Decomp_RCF_input

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Field_Type_Mod, Only : &
    field_type

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_lsm,             &
    stashcode_prog_sec

Use Field_types

Use UM_ParParams

USE UM_ParVars, ONLY : &
    change_decomposition

USE cppxref_mod, ONLY: &
    ppx_atm_compressed

IMPLICIT NONE

! Arguments
Integer, Intent(In)               :: field_count
Type (field_type), Intent(InOut)  :: fields( field_count )
Type (UM_header_type), Intent(In) :: Dump_hdr

! Local Data
Integer                           :: I
Integer                           :: msg
Integer                           :: Pos
Integer                           :: orig_decomp
Integer                           :: ErrorStatus
Character (Len=80)                :: Cmessage
Character (Len=*), Parameter      :: RoutineName = 'Rcf_ReadLSMIn'

! External Routines
External GC_IBcast

!------------------------------------------------------------------
! Set decomposition to be input - communications later require it
!------------------------------------------------------------------
orig_decomp = current_decomp_type
If ( orig_decomp /= decomp_rcf_input ) Then
  Call Change_Decomposition( decomp_rcf_input )
End If

!-------------------------------------------------------------------
! Check that relevant memory is available
!-------------------------------------------------------------------
If ( .NOT. Allocated( Local_LSM_In) ) Then
  Cmessage = 'Local LSM space not allocated!'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If (.NOT. Allocated( Glob_LSM_In ) ) Then
  Cmessage = 'global LSM space not allocated!'
  ErrorStatus = 30
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-------------------------------------------------------------------
! Find Land-Sea Mask amongst input fields
! Allow return code of pos = 0
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields, field_count, POS, .TRUE. )

If (pos == 0) Then
  ! Default - set land sea mask to .FALSE. everywhere
  Local_LSM_In(:) = .FALSE.
  Glob_LSM_In(:)  = .FALSE.

Else
!------------------------------------------------------------------
! Read Local LSM on all PEs
!------------------------------------------------------------------
  ! allocate space in data type and copy data to Local_Lsm_In
  Call Rcf_Alloc_Field( fields(pos) )
  Call Rcf_Read_Field( fields(pos), Dump_hdr, decomp_rcf_input )
  Local_Lsm_In(:) = fields(pos) % Data_Log(:,1)
  Call Rcf_Dealloc_Field( fields(pos) )

!--------------------------------------------------------------------
! And need global LSM
!-------------------------------------------------------------------
  Call Rcf_Gather_Field( local_lsm_in,          glob_lsm_in,        &
                         fields(pos) % row_len, fields(pos) % rows, &
                         fields(pos) % glob_row_len,                &
                         fields(pos) % glob_rows, 0,                &
                         gc_all_proc_group )

! Note that this communication should ideally use GC_BBcast and
! use Kind to determine number of bytes, but we'll not bother.
! After all Gather_Field (above) and the I/O assume that a
! Logical is the same size as a Real as is an Integer. We will do
! likewise.

  msg = 802
  Call GC_IBcast( msg, fields(pos) % glob_level_size, 0, nproc,  &
                  ErrorStatus, glob_lsm_in )

  If ( ErrorStatus /= 0 ) Then
    Cmessage = 'Problem broadcasting global land-sea mask from PE0'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If
!-------------------------------------------------------------------
! Count the local number of land points for addressing
!-------------------------------------------------------------------

local_land_in = 0

Do i = 1, lasize(1,fld_type_p,halo_type_single) * &
    lasize(2,fld_type_p,halo_type_single)
 If ( Local_Lsm_In(i) ) Then
   local_land_in = local_land_in + 1
 End If
End Do

glob_land_in = 0
Do i = 1, glsizep(1) * glsizep(2)
  If ( Glob_Lsm_In(i) ) Then
    glob_land_in = glob_land_in + 1
  End If
End Do

!-----------------------------------------------------------------
! Adjust the land-only fields sizes to match the above sizes
!-----------------------------------------------------------------
Do i = 1, field_count
  If ( fields( i ) % stashmaster % grid_type == ppx_atm_compressed) Then
    If (pos == 0) Then       ! have land only points but no lsm to match
      ErrorStatus = 40
      Cmessage = 'Land-Only fields exist in input dump but *NO* '//&
        'Land-Sea Mask'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    fields( i ) % level_size = local_land_in
  End If
End Do

!-----------------------------------------------------------------
! Reset to original decomposition
!-----------------------------------------------------------------
If ( orig_decomp /= current_decomp_type ) Then
  Call Change_Decomposition( orig_decomp )
End If

Return
End Subroutine Rcf_ReadLSMIn

End Module Rcf_ReadLSMIn_Mod
