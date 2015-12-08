! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reinitialises area, lquid and frozen cloud fractions.

MODULE rcf_derv_cloudfrac_mod

!  Subroutine Rcf_Derv_CloudFrac

! Description:
!   Derive area, liquid and frozen cloud fractions from bulk cf

! Method:
!   A basic copy with no calculation (area) or basic calculation
!   from qcl and qcf (liquid/frozen).

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.4 programming standards.

CONTAINS

SUBROUTINE rcf_derv_cloudfrac( stash_item, fields_out,                         &
                               field_count_out, hdr_out,                       &
                               cloudfrac, data_source )

USE rcf_locate_mod, ONLY :                                                     &
  rcf_locate

USE rcf_alloc_field_mod, ONLY :                                                &
  rcf_alloc_field,                                                             &
  rcf_dealloc_field

USE rcf_read_field_mod, ONLY :                                                 &
  rcf_read_field

USE rcf_data_source_mod, ONLY :                                                &
  data_source_type,                                                            &
  already_processed

USE rcf_stashcodes_mod, ONLY :                                                 &
  stashcode_prog_sec,                                                          &
  stashcode_qcl,       stashcode_qcf,                                          &
  stashcode_bulk_cf,   stashcode_area_cf,                                      &
  stashcode_liquid_cf, stashcode_frozen_cf

USE rcf_field_type_mod, ONLY :                                                 &
  field_type

USE rcf_grid_type_mod, ONLY :                                                  &
  output_grid

USE rcf_level_code_mod, ONLY :                                                 &
  rcf_level_code

USE rcf_umhead_mod, ONLY :                                                     &
  um_header_type

USE printstatus_mod, ONLY :                                                    &
  printstatus,                                                                 &
  prstatus_normal

USE um_parvars, ONLY :                                                         &
  mype

USE decomp_params, ONLY :                                                      &
  decomp_rcf_output

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(in) :: hdr_out
TYPE( field_type ), INTENT(inout), TARGET :: cloudfrac
TYPE( data_source_type ), POINTER :: data_source( : )
INTEGER, INTENT(in)               :: stash_item
INTEGER, INTENT(in)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  area,qcl,qcf

INTEGER                           ::  pos   ! position in array
INTEGER                           ::  i,j,k ! loop index
INTEGER                           ::  start_level, end_level

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  IF ( stash_item == stashcode_bulk_cf ) THEN
    WRITE (6,*) 'Reinitialising Bulk CF as Area CF'
  ELSE IF ( stash_item == stashcode_liquid_cf ) THEN
    WRITE (6,*) 'Reinitialising Liquid CF as fraction of Area CF'
  ELSE IF ( stash_item == stashcode_frozen_cf ) THEN
    WRITE (6,*) 'Reinitialising Frozen CF as fraction of Area CF'
  END IF
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in where available
!----------------------------------------------------------------------
! Bulk Cloud fraction and QCL/QCF if required; will abort if not found

CALL rcf_locate(stashcode_prog_sec, stashcode_area_cf,                         &
                fields_out, field_count_out, pos)
area => fields_out(pos)
CALL rcf_alloc_field( area )
CALL rcf_read_field( area, hdr_out, decomp_rcf_output )

IF ( stash_item == stashcode_liquid_cf .OR.                                    &
     stash_item == stashcode_frozen_cf ) THEN

  CALL rcf_locate(stashcode_prog_sec, stashcode_qcl,                           &
                  fields_out, field_count_out, pos)
  qcl => fields_out(pos)
  CALL rcf_alloc_field( qcl )
  CALL rcf_read_field( qcl, hdr_out, decomp_rcf_output )

  CALL rcf_locate(stashcode_prog_sec, stashcode_qcf,                           &
                  fields_out, field_count_out, pos)
  qcf => fields_out(pos)
  CALL rcf_alloc_field( qcf )
  CALL rcf_read_field( qcf, hdr_out, decomp_rcf_output )
END IF

! ---------------------------------------------------------------------
! Check levels for anything ENDGame related!
! ---------------------------------------------------------------------
start_level = cloudfrac % bottom_level
end_level   = cloudfrac % top_level

! ---------------------------------------------------------------------
! Calculate cloud fractions
! ---------------------------------------------------------------------

! Initialise cloud fraction to be area cloud fraction
! If ENDGame grid, i.e. includes surface level 0, then need to make sure
!  that copy from area cf (no surface level) is done properly.
!  In this case level 0 is copy of level 1.
IF ( start_level == 0 ) THEN
  cloudfrac % data(:,2:end_level+1) = area % data(:,:)
  cloudfrac % data(:,1) = cloudfrac % data(:,2)
ELSE
  cloudfrac % data = area % data
END IF

! If liquid cloud fraction then partition as fraction of qcl
IF ( stash_item == stashcode_liquid_cf ) THEN
  DO k = 1, qcl % levels
    DO i = 1, qcl % level_size
      IF (qcl % data(i,k) > 0.0) THEN
        cloudfrac % data(i,k) = cloudfrac % data(i,k) *                        &
          qcl % data(i,k) / ( qcl % data(i,k) + qcf % data(i,k) )
      ELSE
        cloudfrac % data(i,k) = 0.0
      END IF
    END DO
  END DO
END IF

! If frozen cloud fraction then partition as fraction of qcf
IF ( stash_item == stashcode_frozen_cf ) THEN
  DO k = 1, qcf % levels
    DO i = 1, qcf % level_size
      IF (qcf % data(i,k) > 0.0) THEN
        cloudfrac % data(i,k) = cloudfrac % data(i,k) *                        &
          qcf % data(i,k) / ( qcl % data(i,k) + qcf % data(i,k) )
      ELSE
        cloudfrac % data(i,k) = 0.0
      END IF
    END DO
  END DO
END IF

CALL rcf_locate(stashcode_prog_sec, stash_item,                                &
                fields_out, field_count_out, pos)

data_source( pos ) % source = already_processed

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------

CALL rcf_dealloc_field( area )
IF ( stash_item == stashcode_liquid_cf .OR.                                    &
     stash_item == stashcode_frozen_cf ) THEN
  CALL rcf_dealloc_field( qcl )
  CALL rcf_dealloc_field( qcf )
END IF

RETURN
END SUBROUTINE rcf_derv_cloudfrac
END MODULE rcf_derv_cloudfrac_mod
