! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs transforms after a field has been interpolated.

Module Rcf_Post_Interp_Transform_Mod

!  Subroutine Rcf_Post_Interp_Transform
!
! Description:
!   Wrapper to perform tranformations/processing on a field after
!   it has been interpolated.
!
! Method:
!   Choice of transform/processing is based on stashcode.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE UM_ParVars, Only : &
    mype,                   &
    nproc

Use Rcf_Field_Type_Mod, Only : &
    field_type

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

IMPLICIT NONE

Private field_set_to_min

Contains


Subroutine Rcf_Post_Interp_Transform( output_field, fields_out, &
                                      field_count_out )

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_P_Exner

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_done

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_cca,                 stashcode_cc_lwp,        &
    stashcode_w,                   stashcode_w_adv,         &
    stashcode_exner,               stashcode_q,             &
    stashcode_qcf,                 stashcode_qcl,           &
    stashcode_area_cf,             stashcode_bulk_cf,       &
    stashcode_liquid_cf,           stashcode_frozen_cf,     &
    stashcode_mean_canopyw,        stashcode_can_water_tile,&
    stashcode_qcf2,                stashcode_qrain,         &
    stashcode_qgraup,              stashcode_qc,            &
    stashcode_qT,                  stashcode_prog_sec,      &
    stashcode_proc_phys_sec,       stashcode_tracer_sec,    &
    stashcode_ukca_sec,            stashcode_3d_cca,        &
    stashcode_p

Use Rcf_Recon_Mod, Only : &
    q_min,                &
    w_zero_start,         &
    w_zero_end

USE nlsizes_namelist_mod, ONLY : &
    model_levels

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Submodel_Mod, Only : &
    Atmos_IM

IMPLICIT NONE

! Arguments
Type( field_type ), Intent( InOut ) :: output_field
Type( field_type ), Pointer         :: fields_out(:)
Integer,            Intent(In)      :: field_count_out

! Local variables
Integer                             :: level
Character (Len=*), Parameter        :: RoutineName='Rcf_Post_Interp_Transform'
Character (Len=80)                  :: Cmessage
Integer                             :: ErrorStatus


!---------------------------------------------------------------
! Only do transforms if interpolation is switched on
!---------------------------------------------------------------
If ( output_field % interp == interp_done ) Then

  Select Case( output_field % stashmaster % section )
  
  Case( stashcode_prog_sec )

    ! Which fields do we wish to apply transforms to?
    Select Case( output_field % stashmaster % item )

    Case( stashcode_exner )
      ! convert interpolated P back to Exner
      If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
        Write (6,'(a25)') 'Converting P to exner'
      End If
      ! This would have been changed back to Exner in pre_interp.
      output_field % stashmaster => Rcf_Exppx( Atmos_IM, stashcode_prog_sec, &
                                               stashcode_p)

      Call Rcf_Conv_P_Exner( output_field )

    Case( stashcode_w, stashcode_w_adv )
      ! Zero first level (surface)
      output_field % Data( :,1 ) = 0.0
      If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
        Write (6,*) ' Setting w to zero, level 0'
      End If

      Do level = w_zero_start+1, w_zero_end+1
        output_field % Data( :, level ) = 0.0
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
          Write (6,*) ' Setting w to zero, level ',level-1
        End If
      End Do

      If (w_zero_end /= model_levels) Then
      output_field % Data( :, output_field % levels ) = 0.0
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
          Write (6,*) ' Setting w to zero, level ',model_levels
        End If
      End If

    Case( stashcode_cca )
      Call field_set_to_min( "Convective Cloud Amount", &
                              output_field )
  
      Call field_set_to_max( "Convective Cloud Amount", &
                              output_field )
 
    Case( stashcode_3d_cca )
      Call field_set_to_min( "3D Convective Cloud Amount", &
                              output_field )
 
      Call field_set_to_max( "3D Convective Cloud Amount", &
                              output_field )

    Case( stashcode_cc_lwp )
      Call field_set_to_min( "Liquid Water Path", &
                              output_field )

    Case( stashcode_q )
      Call field_set_to_min( "Specific Humidity", &
                              output_field,       &
                              field_min=q_min )

    Case( stashcode_qcf )
      Call field_set_to_min( "QCF", &
                              output_field )

    Case( stashcode_qcl )
      Call field_set_to_min( "QCL", &
                              output_field )

    Case( stashcode_qcf2)
      Call field_set_to_min( "QCF2", &
                              output_field )

    Case( stashcode_qrain)
      Call field_set_to_min( "QRain", &
                              output_field )

    Case( stashcode_qgraup)
      Call field_set_to_min( "Qgraup", &
                              output_field )

    Case( stashcode_area_cf )
      Call field_set_to_min( "Area Cloud Fraction", &
                              output_field )
 
      Call field_set_to_max( "Area Cloud Fraction", &
                              output_field )

    Case( stashcode_liquid_cf )
      Call field_set_to_min( "Liquid Cloud Fraction", &
                              output_field )
 
      Call field_set_to_max( "Liquid Cloud Fraction", &
                              output_field )

    Case( stashcode_bulk_cf )
      Call field_set_to_min( "Bulk Cloud Fraction", &
                              output_field )
 
      Call field_set_to_max( "Bulk Cloud Fraction", &
                              output_field )

    Case( stashcode_frozen_cf )
      Call field_set_to_min( "Frozen Cloud Fraction", &
                              output_field )
      
      Call field_set_to_max( "Frozen Cloud Fraction", &
                              output_field )

    Case( stashcode_mean_canopyw )
      Call field_set_to_min( "Canopy Water", &
                              output_field )

    Case( stashcode_can_water_tile )
      Call field_set_to_min( "Canopy Water on Tiles", &
                              output_field )

    End Select

  Case( stashcode_proc_phys_sec )

    ! Which fields do we wish to apply transforms to?
    Select Case( output_field % stashmaster % item )

    Case( MOD(stashcode_qc, 1000) )
      Call field_set_to_min( "Cloud water content (qc)", &
                              output_field )

    Case( MOD(stashcode_qT, 1000) )
      Call field_set_to_min( "Total specific humidity (qT)", &
                              output_field )

    End Select

  End Select
End If

Return
End Subroutine Rcf_Post_Interp_Transform


Subroutine field_set_to_min( field_name, output_field, &
                             field_min)
! Iterates across whole of output_field, resetting any
! values below the minimum to be at the minimum.
! Default minimum is zero.

! Arguments
CHARACTER(LEN=*),     Intent( In )      :: field_name
Type( field_type ), Intent( InOut )   :: output_field

!Optional Argument
Real, Optional,     Intent( In )      :: field_min

! Local variables
Integer                               :: i
Integer                               :: j
Integer                               :: count
Integer                               :: istat

Real                                  :: min

! Default minimum is zero
If( Present( field_min ) ) Then
   min = field_min
Else
   min = 0.0
Endif

! Check for fields too low and reset to minimum.
Do j = 1, output_field % levels
   count = 0
   Do i = 1, output_field % level_size
      If ( output_field % Data(i,j) < min ) Then
         output_field % Data(i,j) = min
         count = count + 1
      End If
   End Do

   Call gc_isum (1, nproc, istat, count)
   If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
      write (6,'(a,i4,1x,a,i8)') ' Level ',j, &
           trim(field_name)//' reset to minimum ',count
   End If

End Do

Return
End Subroutine field_set_to_min

Subroutine field_set_to_max( field_name, output_field, &
                             field_max)
! Iterates across whole of output_field, resetting any
! values above the maximum to be at the maximum.
! Default maximum is one.

! Arguments
CHARACTER(LEN=*),     Intent( In )      :: field_name
Type( field_type ), Intent( InOut )   :: output_field

!Optional Argument
Real, Optional,     Intent( In )      :: field_max

! Local variables
Integer                               :: i
Integer                               :: j
Integer                               :: count
Integer                               :: istat

Real                                  :: max

! Default maximum is one
If( Present( field_max ) ) Then
  max = field_max
Else
  max = 1.0
Endif

! Check for fields too high and reset to maximum.
Do j = 1, output_field % levels
  count = 0
  Do i = 1, output_field % level_size
    If ( output_field % Data(i,j) > max ) Then
      output_field % Data(i,j) = max
      count = count + 1
    End If
  End Do

  Call gc_isum (1, nproc, istat, count)
  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    write (6,'(a,i4,1x,a,i8)') ' Level ',j, &
    trim(field_name)//' reset to maximum ',count
  End If

End Do

Return
End Subroutine field_set_to_max

End Module Rcf_Post_Interp_Transform_Mod
