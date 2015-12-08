! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Checks that cloud fraction fields are consistent with q fields.

MODULE Rcf_Cloud_Frac_Chk_Mod

!  Subroutine Rcf_Cloud_Frac_Chk
!
! Description:
!   Ensures that cloud fraction fields are consistent with
!   the humidity (qcf and qcl) fields.
!
! Method:
!   Checks carried out are :-
!   1. Ensure that frozen cloud fraction is zero if QCF is zero.
!   2. Ensure that liquid cloud fraction is zero if QCL is zero.
!   3. Ensure that both area & bulk cloud fractions are zero if both
!      frozen and liquid cloud fractions are zero.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE Rcf_Cloud_Frac_Chk( fields, field_count, grid, decomp, hdr )

USE Rcf_Field_Type_Mod, ONLY : &
    field_type

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY : &
    grid_type

USE Rcf_Locate_Mod, ONLY : &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY : &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY : &
    Rcf_Write_Field

USE Rcf_Alloc_Field_Mod, ONLY : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_Set_Interp_Flags_Mod, ONLY : &
    interp_done

USE Rcf_Stashcodes_Mod, ONLY : &
    stashcode_qcf,             &
    stashcode_qcl,             &
    stashcode_qcf2,            &
    stashcode_area_cf,         &
    stashcode_bulk_cf,         &
    stashcode_frozen_cf,       &
    stashcode_liquid_cf,       &
    stashcode_prog_sec

USE UM_ParVars, ONLY : &
    nproc,                  &
    mype

USE PrintStatus_mod, ONLY : &
    PrintStatus,                &
    PrStatus_Diag

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)
TYPE( grid_type ), INTENT(IN)      :: grid
TYPE( um_header_type ), INTENT(IN) :: hdr
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp

! LOCAL VARIABLES
INTEGER                            :: pos_qcf
INTEGER                            :: pos_qcl
INTEGER                            :: pos_qcf2
INTEGER                            :: pos_cf_a
INTEGER                            :: pos_cf_b
INTEGER                            :: pos_cf_f
INTEGER                            :: pos_cf_l
INTEGER                            :: cldfr_a_changed
INTEGER                            :: cldfr_b_changed
INTEGER                            :: cldfr_f_changed
INTEGER                            :: cldfr_l_changed

INTEGER                            :: i,k
INTEGER                            :: istat
INTEGER                            :: count
INTEGER                            :: start_level
TYPE( field_type ), POINTER        :: qcf
TYPE( field_type ), POINTER        :: qcf2
TYPE( field_type ), POINTER        :: qcl
TYPE( field_type ), POINTER        :: cldfr_a
TYPE( field_type ), POINTER        :: cldfr_b
TYPE( field_type ), POINTER        :: cldfr_f
TYPE( field_type ), POINTER        :: cldfr_l

! Note that any negative qcf, qcl and cloud fractions are reset
! to zero earlier in Rcf_post_process_atmos

!-------------------------------------------------------------------
! Loacte qcl, qcf and cloud fractions in output dump
!-------------------------------------------------------------------

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_qcl,       &
                  fields, field_count, pos_qcl, .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_qcf,       &
                  fields, field_count, pos_qcf, .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_qcf2,      &
                  fields, field_count, pos_qcf2, .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_area_cf,   &
                  fields, field_count, pos_cf_a, .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_bulk_cf,   &
                  fields, field_count, pos_cf_b, .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_frozen_cf, &
                  fields, field_count, pos_cf_f, .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_liquid_cf, &
                  fields, field_count, pos_cf_l, .TRUE. )

!--------------------------------------------------------------------
! Check that frozen cloud fraction is zero if qcf is zero.
!--------------------------------------------------------------------

! Need to check that level size == for qcl and clf_l ?

cldfr_f_changed = 0

IF (pos_qcf /= 0 .AND. pos_cf_f /= 0 ) THEN

  qcf => fields( pos_qcf )
  CALL Rcf_Alloc_Field( qcf )
  CALL Rcf_Read_Field( qcf, hdr, decomp )

  cldfr_f => fields( pos_cf_f )
  CALL Rcf_Alloc_Field( cldfr_f )
  CALL Rcf_Read_Field( cldfr_f, hdr, decomp )

  IF (pos_qcf2 /= 0) THEN
    qcf2 => fields( pos_qcf2 )
    CALL Rcf_Alloc_Field( qcf2 )
    CALL Rcf_Read_Field( qcf2, hdr, decomp )
  END IF

  IF (qcf     % interp == interp_done .OR.   &
      cldfr_f % interp == interp_done ) THEN

    DO k = 1, qcf % levels
      count = 0
      DO i = 1, qcf % level_size

        IF (pos_qcf2 /= 0) THEN ! qcf2 field is present
                                ! so reset cloud fraction only
                                ! if both qcf and qcf2 are zero
          IF ( qcf     % data (i,k) == 0.0  .AND.   &
               qcf2    % data (i,k) == 0.0  .AND.   &
               cldfr_f % data (i,k) >  0.0 ) THEN

            cldfr_f % data (i,k) = 0.0
            cldfr_f_changed = 1
            count = count + 1

          END IF

        ELSE ! qcf2 field not present, only use qcf

          IF ( qcf     % data (i,k) == 0.0  .AND.   &
               cldfr_f % data (i,k) >  0.0 ) THEN

            cldfr_f % data (i,k) = 0.0
            cldfr_f_changed = 1
            count = count + 1

          END IF

        END IF ! on pos_qcf2


      END DO

      CALL gc_isum (1, nproc, istat, count)
      IF (PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
        WRITE (6,'(A,I4,A,I0)') ' Level ',k, &
        ' Cloud fracs (frozen) reset to zero ',count
      END IF

    END DO

  END IF

  IF (pos_qcf2 /= 0 ) THEN
    CALL Rcf_Dealloc_Field( qcf2 )   ! qcf2 no longer required
  END IF
  
  CALL Rcf_Dealloc_Field( qcf )   ! qcf no longer required

END IF


!--------------------------------------------------------------------
! Check that liquid cloud fraction is zero if qcl is zero.
!--------------------------------------------------------------------

! Need to check that level size == for qcl and clf_l ?

cldfr_l_changed = 0

IF (pos_qcl /= 0 .AND. pos_cf_l /= 0 ) THEN

  qcl => fields( pos_qcl )
  CALL Rcf_Alloc_Field( qcl )
  CALL Rcf_Read_Field( qcl, hdr, decomp )

  cldfr_l => fields( pos_cf_l )
  CALL Rcf_Alloc_Field( cldfr_l )
  CALL Rcf_Read_Field( cldfr_l, hdr, decomp )

  IF (qcl     % interp == interp_done .OR.   &
      cldfr_l % interp == interp_done ) THEN

    DO k = 1, qcl % levels
      count = 0
      DO i = 1, qcl % level_size

        IF ( qcl     % data (i,k) == 0.0  .AND.   &
             cldfr_l % data (i,k) >  0.0 ) THEN

          cldfr_l % data (i,k) = 0.0
          cldfr_l_changed = 1
          count = count + 1

        END IF

      END DO

      CALL gc_isum (1, nproc, istat, count)
      IF (PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
        WRITE (6,'(A,I4,A,I0)') ' Level ',k, &
        ' Cloud fracs (liquid) reset to zero ',count
      END IF

    END DO

  END IF
  
  CALL Rcf_Dealloc_Field( qcl )   ! qcl no longer required

END IF


!--------------------------------------------------------------------
! Check that bulk cloud fractions are set to zero
! if both frozen and liquid cloud fractions are zero.
!--------------------------------------------------------------------

cldfr_b_changed = 0

! If the following condition is true we have all the fields we need.
IF (pos_qcf /= 0 .AND. pos_cf_f /= 0 .AND. &
    pos_qcl /= 0 .AND. pos_cf_l /= 0 ) THEN
  IF (pos_cf_b /= 0 ) THEN

    cldfr_b => fields( pos_cf_b )
    CALL Rcf_Alloc_Field( cldfr_b )
    CALL Rcf_Read_Field( cldfr_b, hdr, decomp )

    IF (qcf     % interp == interp_done .OR.   &
        qcl     % interp == interp_done .OR.   &
        cldfr_b % interp == interp_done .OR.   &
        cldfr_f % interp == interp_done .OR.   &
        cldfr_l % interp == interp_done ) THEN

      DO k = 1, qcl % levels
        count = 0
        DO i = 1, qcl % level_size

          IF ( cldfr_f % data (i,k) == 0.0  .AND.   &
               cldfr_l % data (i,k) == 0.0  .AND.   &
               cldfr_b % data (i,k) >  0.0 ) THEN

            cldfr_b % data (i,k) = 0.0
            cldfr_b_changed = 1
            count = count + 1

          END IF

        END DO

        CALL gc_isum (1, nproc, istat, count)
        IF (PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
          WRITE (6,'(A,I4,A,I0)') ' Level ',k, &
          ' Cloud fracs (bulk) reset to zero ',count
        END IF

      END DO

    END IF

  END IF
END IF


!--------------------------------------------------------------------
! Check that area cloud fractions are set to zero
! if both frozen and liquid cloud fractions are zero.
!--------------------------------------------------------------------

start_level = qcl % bottom_level

cldfr_a_changed = 0

! If the following condition is true we have all the fields we need.
IF (pos_qcf /= 0 .AND. pos_cf_f /= 0 .AND. &
    pos_qcl /= 0 .AND. pos_cf_l /= 0 ) THEN
  IF (pos_cf_a /= 0 ) THEN

    cldfr_a => fields( pos_cf_a )
    CALL Rcf_Alloc_Field( cldfr_a )
    CALL Rcf_Read_Field( cldfr_a, hdr, decomp )

    IF (qcf     % interp == interp_done .OR.   &
        qcl     % interp == interp_done .OR.   &
        cldfr_a % interp == interp_done .OR.   &
        cldfr_f % interp == interp_done .OR.   &
        cldfr_l % interp == interp_done ) THEN

      DO k = 1, qcl % levels + start_level - 1
        count = 0
        DO i = 1, qcl % level_size

          IF ( cldfr_f % data (i,k-start_level+1) == 0.0  .AND.   &
               cldfr_l % data (i,k-start_level+1) == 0.0  .AND.   &
               cldfr_a % data (i,k) >  0.0 ) THEN

            cldfr_a % data (i,k) = 0.0
            cldfr_a_changed = 1
            count = count + 1

          END IF

        END DO

        CALL gc_isum (1, nproc, istat, count)
        IF (PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
          write (6,'(A,I4,A,I0)') ' Level ',k, &
          ' Cloud fracs (area) reset to zero ',count
        END IF

      END DO

    END IF

  END IF
END IF

!---------------------------------------------------------------------
! Synchronise `changed' flags
!---------------------------------------------------------------------
CALL GC_Imax( 1, nproc, istat, cldfr_a_changed )
CALL GC_Imax( 1, nproc, istat, cldfr_b_changed )
CALL GC_Imax( 1, nproc, istat, cldfr_f_changed )
CALL GC_Imax( 1, nproc, istat, cldfr_l_changed )

!---------------------------------------------------------------------
! Write out changed fields
!---------------------------------------------------------------------

IF (cldfr_a_changed == 1) THEN
  Call Rcf_Write_Field( cldfr_a, hdr, decomp )
END IF
IF (cldfr_b_changed == 1) THEN
  Call Rcf_Write_Field( cldfr_b, hdr, decomp )
END IF
IF (cldfr_f_changed == 1) THEN
  Call Rcf_Write_Field( cldfr_f, hdr, decomp )
END IF
IF (cldfr_l_changed == 1) THEN
  Call Rcf_Write_Field( cldfr_l, hdr, decomp )
END IF

IF (pos_cf_f /= 0 .AND. pos_cf_l /= 0 ) THEN
  IF (pos_cf_a /= 0 .AND. pos_cf_b /= 0 ) THEN
    CALL Rcf_Dealloc_Field( cldfr_a )
    CALL Rcf_Dealloc_Field( cldfr_b )
  END IF
END IF

IF (pos_qcf /= 0 .AND. pos_cf_f /= 0 ) THEN
  Call Rcf_Dealloc_Field( cldfr_f )
END IF

IF (pos_qcl /= 0 .AND. pos_cf_l /= 0 ) THEN
  Call Rcf_Dealloc_Field( cldfr_l )
END IF

RETURN
END SUBROUTINE Rcf_Cloud_Frac_Chk
END MODULE Rcf_Cloud_Frac_Chk_Mod
