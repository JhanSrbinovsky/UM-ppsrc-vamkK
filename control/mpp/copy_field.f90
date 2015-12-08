! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine COPY_FIELD
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

SUBROUTINE copy_field(                                            &
         field_in, field_out                                      &
,        row_length_in, row_length_out, rows_in, rows_out         &
,        levels_in, levels_out, level_start, level_end            &
,        haloi_in, haloj_in, haloi_out, haloj_out                 &
,        fld_type, l_haloes, l_swap, l_vector)

! Purpose:
!     This routine copies one field into another, allowing for a
!     different halo size in the two fields.
!    If L_haloes is true then haloes are copied explicitly BUT NB
!    NB input haloes must be greater or equal to output haloes
!    If L_haloes is false then only non-halo regions copied
!    If L_SWAP is true the halos on the destination field will be
!    updated using the standard swap_bounds

! Method:
!          Is described in ;


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


USE ereport_mod, ONLY: ereport
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) :: row_length_in
INTEGER, INTENT(IN) :: row_length_out
INTEGER, INTENT(IN) :: rows_in
INTEGER, INTENT(IN) :: rows_out  
INTEGER, INTENT(IN) :: levels_in
INTEGER, INTENT(IN) :: levels_out
INTEGER, INTENT(IN) :: haloi_in
INTEGER, INTENT(IN) :: haloj_in
INTEGER, INTENT(IN) :: haloi_out
INTEGER, INTENT(IN) :: haloj_out
INTEGER, INTENT(IN) :: level_start
INTEGER, INTENT(IN) :: level_end 
INTEGER, INTENT(IN) :: fld_type     ! Data type for swap bounds

LOGICAL, INTENT(IN) :: l_haloes     ! If true fill haloes explicitly
LOGICAL, INTENT(IN) :: l_swap       ! If true use swap bounds to fill haloes
LOGICAL, INTENT(IN) :: l_vector     ! Vector switch for swap bounds

! Primary Arrays
REAL, INTENT(IN) :: field_in(1-haloi_in: row_length_in+haloi_in,       &
                             1-haloj_in: rows_in+haloj_in, levels_in) 

REAL, INTENT(OUT) :: field_out(1-haloi_out: row_length_out+haloi_out,  &
                         1-haloj_out: rows_out+haloj_out, levels_out)

! Local Variables.
INTEGER  ::  i, j, k  ! Loop indices
INTEGER  ::  levels

INTEGER  ::  ErrorStatus ! ErrorStatus for ereport

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'copy_field'
CHARACTER(LEN=80)             :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('COPY_FIELD',zhook_in,zhook_handle)
IF (l_haloes) THEN
! ----------------------------------------------------------------------
! Section 1.  If L_haloes = .true. then explicitly copy halo region
!             In this case halo_in must be >= halo_out
! ----------------------------------------------------------------------

  IF (       haloi_out  >   haloi_in                              &
       .OR.  haloj_out  >   haloj_in) THEN
    cmessage = "Error in COPY_FIELD: Output halo must be smaller "&
             //"than input halo"
    ErrorStatus = 1
    CALL ereport(RoutineName, ErrorStatus, cmessage)
  END IF

  DO k = level_start, level_end
    DO j = 1-haloj_out, rows_out+haloj_out
      DO i = 1-haloi_out, row_length_out+haloi_out
        field_out(i,j,k) = field_in(i,j,k)
      END DO
    END DO
  END DO

ELSE
! ----------------------------------------------------------------------
! Section 2.  If L_ haloes= .false. then only copy non-halo region
!              To fill haloes set L_swap = .true to use swap bounds
! ----------------------------------------------------------------------
  DO k = level_start, level_end
    DO j = 1, rows_out
      DO i = 1, row_length_out
        field_out(i,j,k) = field_in(i,j,k)
      END DO
    END DO
  END DO

  IF (l_swap) THEN

    levels = level_end - level_start + 1

! DEPENDS ON: swap_bounds
    CALL swap_bounds(                                             &
                   field_out(1-haloi_out,1-haloj_out,level_start) &
,                  row_length_out, rows_out, levels               &
,                  haloi_out, haloj_out, fld_type, l_vector)
  END IF !  L_swap

END IF ! L_haloes

IF (lhook) CALL dr_hook('COPY_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE copy_field

