! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads in auxillary data into dump

Module Rcf_Aux_File_Mod

!  Subroutine Rcf_Aux_File  - reads auxillary data
!
! Description:
!   Reads data from external dumps and incorporates it into the
!   output dump. Four modes are supported
!                Tracer Data
!                User Prognostics
!                UARS Data
!                Area Tranplants
!
! Method:
!   A fields array is set up for the auxillary file and the relevant
!   fields located therein. Data is copied as appropriate for the
!   Mode (ie transplants within an area, all of a user prognostic
!   copied, upper levels only of UARS or Tracers copied).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! Parameters describing actions possible
Integer, Parameter       :: tracers   = 1
Integer, Parameter       :: user_prog = 2
Integer, Parameter       :: uars_data = 3
Integer, Parameter       :: transplant= 4

Contains

Subroutine Rcf_Aux_File( Hdr_Aux, Hdr_Out, Fields_Out, Field_Count_Out,&
                         Action, Sctn_code_aux, Item_Code_aux,         & 
                         Sctn_code_out, Item_code_out )

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

USE mask_compression, ONLY: compress_to_mask

Use Rcf_UMhead_Mod, Only : &
    Um_Header_Type

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

USE UM_ParVars, Only : &
    mype

Use rcf_trans_mod              ! all of it

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Setup_Field_mod, Only : &
    Rcf_Setup_Field

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_ReadUMhdr_Mod, Only : &
    Rcf_ReadUMhdr

Use Rcf_HeadAddress_Mod, Only :  &
IC_XLen,                         IC_YLen,    &
FH_DTYear,                       FH_VTYear,  &
FH_DTDayNo,                      FH_VTDayNo,       &
FH_Dataset,                      FH_Dataset_Ancil

Use Rcf_Grid_Type_Mod, Only :&
    Output_Grid

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only : &
    local_land_out,     &
    local_lsm_out

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Global_To_Local_Subdomain

Use cppxref_mod, Only :             &
    ppx_type_real,                  &
    ppx_type_int,                   &
    ppx_type_log,                   &
    ppx_atm_compressed

IMPLICIT NONE

! Arguments
Type( Um_Header_type ), Intent(InOut)  :: Hdr_Aux ! only unit no.,
                                                  ! file already open
Type( Um_Header_type ), Intent(In)     :: Hdr_Out
Type( Field_Type ), Pointer            :: Fields_Out(:)

Integer, Intent(In)                    :: Field_Count_Out
Integer, Intent(In)                    :: sctn_code_aux
Integer, Intent(In)                    :: item_code_Aux
Integer, Intent(In)                    :: sctn_code_out
Integer, Intent(In)                    :: item_code_out ! only for UP
Integer, Intent(In)                    :: Action

! Local variables
Integer                     :: i, j    ! loopers
Integer                     :: k,l     ! loopers
Integer                     :: pos_out ! field position
Integer                     :: pos_aux ! field position
Integer                     :: lrow1   ! Local positions of
Integer                     :: lrow2   ! trans data validity
Integer                     :: lcol1
Integer                     :: lcol2
Integer                     :: level_base
Integer                     :: level_top
Integer                     :: size
Integer                     :: field_count_aux
Integer                     :: ErrorStatus
Integer                     :: Copy_Count ! counter for no. of x a field
                                          ! is copied

Integer, Parameter          :: st_no_data = -3

Character (Len=20)          :: title
Character (Len=*), Parameter:: RoutineName='Rcf_Aux_File'
Character (Len=80)          :: Cmessage

Type( field_type ), Pointer :: fields_aux(:)

Nullify( fields_aux )
!------------------------------------------------------------------
! Read in the Auxillary file header
!------------------------------------------------------------------
Call Rcf_ReadUMhdr( Hdr_Aux )

!------------------------------------------------------------------
! If uars or tracers check that data time of the aux file matches
! verification time of the output file
!------------------------------------------------------------------
If ( action == tracers .OR. action == uars_data ) Then
  Do i = 0, 6
    If ( Hdr_Aux % FixHd( FH_DTYear + i ) /= &
         Hdr_Out % FixHd( FH_VTYear + i ) ) Then
      Write (6,*) 'Mismatch in date info for auxillary file'
      Write (6,*) 'Aux file Data times = ', (Hdr_Aux % Fixhd( j ), &
                                             j = FH_DTYear, FH_DTDayNo )
      Write (6,*) 'Out file Ver. times = ', (Hdr_Out % Fixhd( j ), &
                                             j = FH_VTYear, FH_VTDayNo )

      Cmessage = 'Date information mismatch between aux and dump files'
      ErrorStatus = 10
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End Do
End If

!-------------------------------------------------------------------
! Check that resolutions of two files match
!-------------------------------------------------------------------
If ( Hdr_Aux % IntC( IC_XLen ) /= Hdr_Out % IntC( IC_XLen ) .OR. &
     Hdr_Aux % IntC( IC_YLen ) /= Hdr_Out % IntC( IC_YLen ) ) Then

  Cmessage = 'Dimensions of AUX file and dump file do not match'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-------------------------------------------------------------------
! Set up field data-types for aux file
! (Grid resolutions should be as for output grid)
!-------------------------------------------------------------------
title = 'Auxillary File'
If ( action == user_prog .AND.                                 &
     Hdr_Aux % FixHd(FH_Dataset)  == FH_Dataset_Ancil ) Then
  ! User prognostic ancillary files have full fields for land-only
  ! fields
  Call Rcf_Setup_Field( fields_aux, Hdr_Aux, Output_Grid,      &
                        field_count_aux, title,                &
                        Output_Grid % loc_p_rows *             &
                        Output_Grid % loc_p_row_length )
Else
  Call Rcf_Setup_Field( fields_aux, Hdr_Aux, Output_Grid,      &
                        field_count_aux, title, local_land_out )
End If

!-------------------------------------------------------------------
! Main data handling/replacement - start with TRANS data
!-------------------------------------------------------------------
If ( action == transplant ) Then
  Do i = 1, num_trans

    If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
      Write (6,*) 'Transplanting data for stashcode ', itemc_array( i )
    End If

    ! Read aux field
    Call Rcf_Locate( sctnc_array( i ), itemc_array( i ),             &
                     fields_aux, field_count_aux, pos_aux)
    Call Rcf_Alloc_Field( fields_aux( pos_aux ) )
    Call Rcf_Read_Field( fields_aux( pos_aux ), Hdr_Aux,             &
                         decomp_rcf_output )

    ! Read dump field
    Call Rcf_Locate( sctnc_array( i ), itemc_array( i ),             & 
                     fields_out, field_count_out, pos_out)

    Call Rcf_Alloc_Field( fields_out( pos_out ) )
    Call Rcf_Read_Field( fields_out( pos_out ), Hdr_Out,             &
                         decomp_rcf_output )

    ! Cannot (yet) do a transplant for a land compressed field
    If (fields_out(pos_out) % stashmaster % grid_type ==             &
                                            ppx_atm_compressed) Then
      Cmessage = 'Cannot transplant a land compressed field'
      ErrorStatus = 30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    ! Convert the co-ordinates for replacement with local ones
    Call Rcf_Global_to_Local_Subdomain( .FALSE., .FALSE.,             &
                    fields_out( pos_out ) % stashmaster % grid_type , &
                    mype,                                             &
                    Row2_array( i ), Col2_array( i ),                 &
                    Row1_array( i ), Col1_array( i ),                 &
                    lrow2, lcol2, lrow1, lcol1 )

    ! If the local area is on my pe, do the transplant
    If ( lrow1 /= st_no_data .AND. lrow2 /= st_no_data .AND. &
         lcol1 /= st_no_data .AND. lcol2 /= st_no_data ) Then

      Select Case ( fields_out( pos_out ) % stashmaster % data_type )
        Case (ppx_type_real)
          Do j = lrow1, lrow2
            Do k = lcol1, lcol2
              l = (j-1) * fields_out(pos_out) % row_len + k
              fields_out(pos_out) %                      &
                    data(l,lev1_array(i):lev2_array(i))= &
              fields_aux(pos_aux) % data(l,lev1_array(i):lev2_array(i))
            End Do
          End Do

        Case (ppx_type_int)
          Do j = lrow1, lrow2
            Do k = lcol1, lcol2
              l = (j-1) * fields_out(pos_out) % row_len + k
              fields_out(pos_out) %                          &
                    data_int(l,lev1_array(i):lev2_array(i))= &
              fields_aux(pos_aux) %                          &
                    data_int(l,lev1_array(i):lev2_array(i))
            End Do
          End Do

        Case (ppx_type_log)
          Do j = lrow1, lrow2
            Do k = lcol1, lcol2
              l = (j-1) * fields_out(pos_out) % row_len + k
              fields_out(pos_out) %                          &
                    data_log(l,lev1_array(i):lev2_array(i))= &
              fields_aux(pos_aux) %                          &
                    data_log(l,lev1_array(i):lev2_array(i))
            End Do
          End Do

      End Select
    End If

    ! Write out the field
    Call Rcf_Write_Field( fields_out( pos_out ), Hdr_Out, &
                         decomp_rcf_output)

    Call Rcf_Dealloc_Field( fields_out( pos_out ) )
    Call Rcf_Dealloc_Field( fields_aux( pos_aux ) )

  End Do

Else       ! not transplant....
!----------------------------------------------------------------
! Loop through the auxillary fields for consideration for
! output dump inclusion.
! After problems with an ocean dump containing a time series diagnostic
! feild with identical section and item no. but different size the
! Copy_Count now forces only the first instance in the dump to be used
!----------------------------------------------------------------
  Copy_Count=0
  Do i = 1, field_count_aux
    If ( (action == tracers     .OR.   &
          action == uars_data   .OR.   &
          action == user_prog ) .AND.  &
         ( item_code_aux == fields_aux( i ) % stashmaster % item .AND. &
          sctn_code_aux == fields_aux( i ) % stashmaster % section ) ) Then

      ! Increment the copy count and check its the 1st instance.
      Copy_Count=Copy_Count + 1
      If ( Copy_Count == 1 ) Then

        If ( action == user_prog ) Then
          Call Rcf_Locate( sctn_code_out, item_code_out,               &
                           fields_out, field_count_out, pos_out)
        Else
          Call Rcf_Locate( fields_aux( i ) % stashmaster % section,    &
                           fields_aux( i ) % stashmaster % item,       &
                           fields_out, field_count_out, pos_out )
        End If

        Call Rcf_Alloc_Field( fields_aux( i ) )
        Call Rcf_Read_Field( fields_aux(i), Hdr_Aux, decomp_rcf_output )

        Call Rcf_Alloc_Field( fields_out( pos_out ) )

!---------------------------------------------------------------
! Copy fields if user prog
!---------------------------------------------------------------
        If ( action == user_prog ) Then

          ! If a land only field from ancillary, can compress onto
          ! output field. This assumes a real field only
          If (fields_aux( i ) % stashmaster % grid_type ==             &
                                ppx_atm_compressed .AND.               &
              Hdr_Aux % FixHd( FH_Dataset) == FH_Dataset_Ancil ) Then

            Do j = 1, fields_out( pos_out ) % levels
              Call compress_to_mask( fields_aux( i ) % Data(:,j),        &
                  fields_out( pos_out ) % Data(:,j),  &
                  local_lsm_out,                      &
                  fields_aux( i ) % level_size,       &
                  size )
            End Do

          Else    ! not land compressed

          Select Case ( fields_out( pos_out ) % stashmaster % data_type)
            Case (ppx_type_real)
              fields_out( pos_out ) % Data( :, : ) =                   &
                                      fields_aux( i ) % Data( :, : )

            Case (ppx_type_int)
              fields_out( pos_out ) % Data_int( :, : ) =               &
                                     fields_aux( i ) % Data_int( :, : )

            Case (ppx_type_log)
              fields_out( pos_out ) % Data_Log( :, : ) =               &
                                     fields_aux( i ) % Data_Log( :, : )

          End Select
          End If

        Else             ! Must be uars or tracers

          ! If levels don't match for tracers, issue a warning
          If ( action == tracers .AND. fields_out( pos_out ) % levels/=&
                                       fields_aux( i ) % levels) Then
            ErrorStatus = -40
            Cmessage = 'Not all tracer levels have been initialised!'
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          Call Rcf_Read_Field( fields_out( pos_out ), Hdr_Out,         &
                                                     decomp_rcf_output)

        ! Only copy the aux levels over the top of the output levels
          level_top  = fields_out( pos_out ) % levels
          level_base = fields_out( pos_out ) % levels -                &
                       fields_aux( i ) % levels + 1

          Select Case ( fields_out( pos_out ) % stashmaster % data_type)
            Case (ppx_type_real)
            fields_out( pos_out ) % Data( :, level_base : level_top) = &
                                      fields_aux( i ) % Data( :, : )

            Case (ppx_type_int)
              fields_out(pos_out) % Data_int( :,level_base:level_top)= &
                                    fields_aux( i ) % Data_int( :, : )

            Case (ppx_type_log)
              fields_out(pos_out) % Data_Log( :,level_base:level_top)= &
                                    fields_aux( i ) % Data_Log( :, : )

          End Select

        End If

        Call Rcf_Write_Field( fields_out( pos_out ), Hdr_Out,          &
                              decomp_rcf_output )

        Call Rcf_Dealloc_Field( fields_out( pos_out ) )
        Call Rcf_Dealloc_Field( fields_aux( i ) )

      Else  ! Copy_Count /= 1, must have already copied field over

        ErrorStatus = -50
        Cmessage="Was about to overwrite user_prog data with 2nd field."
        Call EReport ( RoutineName, ErrorStatus, Cmessage )

      End If ! Copy_Count = 1
    End If
  End Do

End If

! Clean up the fields
Deallocate( fields_aux )
Nullify( fields_aux )

Return
End Subroutine Rcf_Aux_File
End Module Rcf_Aux_File_Mod
