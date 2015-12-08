! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the orography for the reconfiguration

MODULE Rcf_Set_Orography_Mod
IMPLICIT NONE

!  Subroutine Rcf_Set_Orography - sets the orographies
!
! Description:
! This module sets up input, output and interpolated orography
! fields for interpolation (if required ) (particularly for height
! field generation.
!
! Method:
!  Read input and output orographies - if required interpolate the
!  input orography. Do orogrphic blending for LAM and check if
!  orographic changes force vertical interpolation.
!
!  Orographic blending blends from "outside to inside" of the domain
!  with the specified weight. The blending zone *includes* the
!  Rim - thus a weight of 1 should be specified to leave the Rim
!  untouched.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!

CONTAINS

SUBROUTINE Rcf_Set_Orography( fields_in, fields_out, field_count_in, &
                              field_count_out, hdr_in, hdr_out,      &
                              data_source, orog_in, orog_out,        &
                              interp_orog )

USE Rcf_Gather_Field_Mod, ONLY : &
    Rcf_Gather_Field_Real

USE Rcf_Scatter_Field_Mod, ONLY : &
    Rcf_Scatter_Field_Real

USE PrintStatus_mod, ONLY : &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParVars, ONLY :      &
    mype,                   &
    nproc,                  &
    datastart,              &
    current_decomp_type,    &
    gc_all_proc_group,      &
    change_decomposition

USE Rcf_Locate_Mod, ONLY : &
    Rcf_Locate

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY : &
    Input_Grid,               &
    Output_Grid

USE Rcf_Field_Type_Mod, ONLY : &
    field_type

USE decomp_params, ONLY : &
    decomp_rcf_input,        &
    decomp_rcf_output

USE Rcf_Interp_Weights_Mod, ONLY : &
    h_int_active

USE Rcf_V_Int_Ctl_Mod, ONLY : &
    v_int_active

USE Rcf_field_equals_mod, ONLY  : &
    Rcf_field_equals

USE Rcf_Set_Interp_Flags_Mod, ONLY :  &
    interp_copy,                      &
    interp_h_only

USE Rcf_Data_Source_Mod, ONLY : &
    data_source_type,           &
    Ancillary_File,             &
    already_processed

USE Rcf_Read_Field_Mod, ONLY : &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY : &
    Rcf_Write_Field

USE Rcf_Interpolate_Mod, ONLY : &
    Rcf_Interpolate

USE Rcf_Alloc_Field_Mod, ONLY : &
    Rcf_Alloc_Field

USE Rcf_Stashcodes_Mod, ONLY : &
    stashcode_orog,            &
    stashcode_prog_sec

USE Rcf_ReadNL_Horizont_Mod, ONLY : &
    orog_blend_width,               &
    blend_weights

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER             :: fields_in(:)
TYPE( field_type ), POINTER             :: fields_out(:)
TYPE( field_type ), POINTER             :: orog_in
TYPE( field_type ), POINTER             :: orog_out
TYPE( field_type ), TARGET              :: interp_orog
TYPE( data_source_type ), INTENT(INOUT) :: data_source(:)
TYPE( um_header_type ), INTENT(IN)      :: hdr_in
TYPE( um_header_type ), INTENT(IN)      :: hdr_out
INTEGER, INTENT(IN)                     :: field_count_in
INTEGER, INTENT(IN)                     :: field_count_out

! Local variables
INTEGER                              :: i
INTEGER                              :: j
INTEGER                              :: k
INTEGER                              :: pos_in
INTEGER                              :: pos_out
INTEGER                              :: decomp_old   ! old decomposition
INTEGER                              :: stat         ! gcom status
INTEGER                              :: v_on = 0     ! vert interp flag
LOGICAL                              :: ancil_orog
REAL                                 :: weight       ! blending weight
REAL, ALLOCATABLE                    :: orog_out_fullfield   ( :, : )
REAL, ALLOCATABLE                    :: orog_interp_fullfield( :, : )
TYPE (field_type)                    :: dummy

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE (6,'(A)') 'Processing Orography (stashcode 33) '
END IF

!------------------------------------------------------------------
! Find output orography
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos_out, .TRUE. )

! Check we have found it.
IF (pos_out == 0) THEN
  WRITE (6,'(A)') 'Output orography not found.  Attempting to continue.'
  ALLOCATE(orog_in)
  NULLIFY(orog_in % data)
  ALLOCATE(orog_out)
  NULLIFY(orog_out % data)
  GOTO 9999
END IF

!---------------------------------------------------------------
! Find origin of orography
!----------------------------------------------------------------

IF ( data_source( pos_out ) % source == Ancillary_File ) THEN
  ancil_orog = .TRUE.
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
    WRITE (6,'(A)') 'Using Ancillary Orography'
  END IF
ELSE
  ancil_orog = .FALSE.
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
    IF (h_int_active) THEN
      WRITE (6,'(A)') 'Using interpolated Orography'
    ELSE
      WRITE (6,'(A)') 'Copying input Orography'
    END IF
  END IF
END IF

!------------------------------------------------------------------
! find and read input orography
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_in, field_count_in, pos_in )
orog_in  => fields_in(pos_in)
CALL Rcf_Alloc_Field( orog_in )
CALL Rcf_Read_Field( orog_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Setup output orography
!------------------------------------------------------------------
orog_out => fields_out(pos_out)
CALL Rcf_Alloc_Field( orog_out )

! Set the sizes of interp_orog to be those of orog_out
CALL rcf_field_equals(interp_orog, orog_out)

!------------------------------------------------------------------
! Setup inteperpolated orography
!------------------------------------------------------------------
IF (h_int_active) THEN
  orog_in % interp = interp_h_only
ELSE
  orog_in % interp = interp_copy
END IF

CALL Rcf_Alloc_Field( interp_orog )
CALL rcf_interpolate( orog_in, interp_orog, Input_Grid,       &
                      Output_Grid, dummy, dummy )

!-------------------------------------------------------------------
! Check for which output orography is required
!-------------------------------------------------------------------
IF (ancil_orog) THEN
  ! read the ancillary back in from output dump
  CALL Rcf_Read_Field( orog_out, Hdr_Out, decomp_rcf_output )

  !----------------------------------------------------------------
  ! Perform the Topog masking
  !----------------------------------------------------------------
  decomp_old = decomp_rcf_output
  IF (current_decomp_type /= decomp_rcf_output) THEN
    decomp_old = current_decomp_type
    CALL Change_Decomposition( decomp_rcf_output )
  END IF

  IF ( orog_blend_width > 0 ) THEN

    ALLOCATE( orog_out_fullfield(    orog_out % glob_row_len, &
                                     orog_out % glob_rows ) )
    ALLOCATE( orog_interp_fullfield( orog_out % glob_row_len, &
                                     orog_out % glob_rows ) )

    ! Gather orogrophies on PE 0
    ! Cannot use generic routine as fullfield is 2D so doesn't match
    ! in rank!
    CALL Rcf_Gather_Field_Real( orog_out % DATA(:,1),                &
                                orog_out_fullfield,                  &
                                orog_out % row_len, orog_out % rows, &
                                orog_out % glob_row_len,             &
                                orog_out % glob_rows, 0,             &
                                gc_all_proc_group )

    CALL Rcf_Gather_Field_Real( interp_orog % DATA(:,1),             &
                                orog_interp_fullfield,               &
                                orog_out % row_len, orog_out % rows, &
                                orog_out % glob_row_len,             &
                                orog_out % glob_rows, 0,             &
                                gc_all_proc_group )

    ! Do the orography blending on PE 0
    IF (mype == 0) THEN

      ! Northern and Southern Strips (including corners)
      DO i = 1, orog_out % glob_row_len
        DO j = 1, orog_blend_width

          ! First determine which weight to use
          ! Western corners
          IF ( i < orog_blend_width ) THEN
            weight = blend_weights( MIN(i,j) )

          ! Eastern corners
          ELSE IF ( i > orog_out % glob_row_len - orog_blend_width + 1)&
                                                                   THEN
            weight = blend_weights(                                 &
                               MIN( orog_out % glob_row_len - i + 1, j))

          ! Middle section
          ELSE
            weight = blend_weights( j )

          END IF

          ! Set the blended field for the Southern strip
          k = j
          orog_out_fullfield(i,k) =                                   &
                    orog_interp_fullfield(i,k) * weight +             &
                    orog_out_fullfield(i,k) * (1.0 - weight )

          ! Set the blended field for the Northern strip
          k = orog_out % glob_rows - j + 1
          orog_out_fullfield(i,k) =                                   &
                    orog_interp_fullfield(i,k) * weight +             &
                    orog_out_fullfield(i,k) * (1.0 - weight )
        END DO
      END DO

      ! Western and Eastern Strips (excluding corners)
      DO i = 1, orog_blend_width
        DO j = orog_blend_width + 1, orog_out % glob_rows -           &
                                     orog_blend_width

          ! Set the weight used
          weight = blend_weights( i )

          ! Set the blended field for the Western Strip
          k = i
          orog_out_fullfield(k,j) =                                   &
                    orog_interp_fullfield(k,j) * weight +             &
                    orog_out_fullfield(k,j) * (1.0 - weight )

          ! Set the blended field for the Eastern Strip
          k = orog_out % glob_row_len - i + 1
          orog_out_fullfield(k,j) =                                   &
                    orog_interp_fullfield(k,j) * weight +             &
                    orog_out_fullfield(k,j) * (1.0 - weight )

        END DO
      END DO

    END IF

    ! Only need to scatter the orog_out_fullfield
    CALL Rcf_Scatter_Field_Real( orog_out % DATA, orog_out_fullfield, &
                                 orog_out % row_len, orog_out % rows, &
                                 orog_out % glob_row_len,             &
                                 orog_out % glob_rows, 0,             &
                                 gc_all_proc_group )

    CALL Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )

    DEALLOCATE( orog_out_fullfield )
    DEALLOCATE( orog_interp_fullfield )

  END IF

  ! Change decomposition back
  IF ( current_decomp_type /= decomp_old ) THEN
    CALL Change_Decomposition( decomp_old )
  END IF

ELSE
  ! Use interpolated orography for output
  orog_out => interp_orog
  CALL Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )
  ! We have now handled the orography.
  data_source( pos_out ) % source = already_processed
END IF

!------------------------------------------------------------------
! If there is a difference between input and output orographies,
! we need to turn on vertical interpolation
!------------------------------------------------------------------
IF ( .NOT. v_int_active ) THEN
  IF ( orog_in % glob_level_size /= orog_out % glob_level_size ) THEN
    IF (ancil_orog) THEN    ! Only switch on interpolation if ancillary
      v_int_active = .TRUE.
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
        WRITE (6,'(A)') 'Vertical interpolation has been switched on '//&
        'due to a change in orography'
      END IF
    END IF
  ELSE
    DO i = 1, orog_in % level_size
      IF ( ABS( orog_in % DATA( i, 1) - orog_out % DATA( i, 1 ) ) > &
           EPSILON( 1.0 ) ) THEN
        v_on = 1
        EXIT
      END IF
    END DO

    ! Need to make sure all PEs turn on v interp if 1 does
    CALL GC_IMAX( 1, nproc, stat, v_on )
    IF ( v_on == 1 ) THEN
      v_int_active = .TRUE.
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
        WRITE (6,'(A)') 'Vertical interpolation has been switched on '//&
        'due to a change in orography'
      END IF
    END IF
  END IF
END IF

9999 CONTINUE
RETURN
END SUBROUTINE Rcf_Set_Orography
END MODULE Rcf_Set_Orography_Mod
