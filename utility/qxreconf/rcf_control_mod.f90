! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ The main control routine for the reconfiguration

Module Rcf_Control_Mod

!  Subroutine Rcf_Control - main control routine
!
! Description:
!    Top level control of the reconfiguration
!
! Method:
!    Headers, Fields land-sea masks and interpolation weights are
!    set up. Then ancillaries are read and main dump creation done.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Control( hdr_in, hdr_out )

Use Rcf_Setup_Header_Mod, Only : &
    Rcf_Setup_Header

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use  Rcf_Setup_Field_mod, Only : &
     Rcf_Setup_Field

Use Rcf_init_h_interp_mod, Only : &
    Rcf_init_h_interp

Use Rcf_ReadLSMIn_mod, Only : &
    Rcf_ReadLSMIn

Use Rcf_Setup_LSM_out_mod, Only : &
    Rcf_setup_lsm_out

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,           &
    Output_Grid

Use Rcf_Ancil_mod, Only : &
    Rcf_Ancil

Use Rcf_Create_Dump_Mod, Only : &
    Rcf_Create_Dump

Use Rcf_WriteUMhdr_Mod, Only : &
    Rcf_WriteUMhdr

USE rcf_set_data_source_mod, ONLY: &
    rcf_set_data_source

USE rcf_data_source_mod, ONLY: &
    data_source_type

Implicit None

! Arguments
Type (um_header_type), Intent(InOut)   :: hdr_in
Type (um_header_type), Intent(InOut)   :: hdr_out

! Local variables
Type (field_type), Pointer :: fields_in(:)
Type (field_type), Pointer :: fields_out(:)
Integer                    :: field_count_in
Integer                    :: field_count_out
Character (Len=20)         :: title
Type (data_source_type), Pointer :: data_source(:)

!---------------------------------------------------------------
! Initialise fields data-type(s)
!---------------------------------------------------------------
Nullify( fields_in  )
Nullify( fields_out )

!----------------------------------------------------------------
! Setup the output header information (including lookups)
!----------------------------------------------------------------
Call Rcf_Setup_Header( hdr_in, hdr_out )

!---------------------------------------------------------------
! Write out the header
!---------------------------------------------------------------
Call Rcf_WriteUMhdr( Hdr_Out )

!----------------------------------------------------------------
! Setup the field data-types for both the input and output grids
!----------------------------------------------------------------
title = 'Input grid'
Call Rcf_Setup_Field( fields_in, Hdr_In, Input_Grid,                &
                      field_count_in, title)

title = 'Output grid'
Call Rcf_Setup_Field( fields_out, Hdr_Out, Output_Grid,             &
                      field_count_out, title )

!------------------------------------------------------------------
! Setup horizontal interpolation weights
!------------------------------------------------------------------
Call Rcf_init_h_interp( Input_Grid, Output_Grid, Hdr_In, Hdr_Out )

!-------------------------------------------------------------------
! Setup input and output land-sea masks and coastal adjustment
! indexes etc.
!-------------------------------------------------------------------
Call Rcf_ReadLSMIn( Hdr_In, fields_in, field_count_in )
Call Rcf_Setup_LSM_Out( hdr_out, fields_in, field_count_in, fields_out, &
                        field_count_out )

!-------------------------------------------------------------------
! Setup data source
!-------------------------------------------------------------------
! set the source of the data for output dump
NULLIFY(data_source)
Call Rcf_Set_Data_Source( data_source, fields_in, fields_out,  &
                          field_count_in, field_count_out,     &
                          hdr_in, hdr_out )

!-------------------------------------------------------------------
! Do the ancillary processing
!-------------------------------------------------------------------
Call Rcf_Ancil ( Hdr_In, Hdr_Out,              &
                 Fields_In, Field_Count_In,    &
                 Fields_Out, Field_Count_Out,  &
                 data_source )

!-------------------------------------------------------------------
! Create the output dump
!-------------------------------------------------------------------
Call Rcf_Create_Dump( hdr_in, hdr_out, fields_in, fields_out,       &
                      field_count_in, field_count_out, data_source )


Return
End Subroutine Rcf_Control
End Module Rcf_Control_Mod
