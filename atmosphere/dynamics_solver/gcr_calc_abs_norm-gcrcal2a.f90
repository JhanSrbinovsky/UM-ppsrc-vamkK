! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Calculate max value of error over all processors
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver


      SUBROUTINE GCR_calc_abs_norm(                                     &
     &                            Error, off_x, off_y,                  &
     &                            n_proc, model_domain,                 &
     &                            at_extremity,                         &
     &                            row_length, rows,                     &
     &                            model_levels, Abs_Norm)



      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, n_proc                                                          &
     &, off_x                                                           &
     &, off_y                                                           &
     &, model_domain

      LOGICAL                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      REAL                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Arguments with Intent OUT. ie: variables Output only
      REAL                                                              &
     &  Abs_norm

! Local variables

      INTEGER i, j, k, info                                             &
     &, i_start, i_end                                                  &
     &, j_start, j_end

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External Routines:

      EXTERNAL                                                          &
     &  gc_rmax

! ----------------------------------------------------------------------
! Section 1.   Calculate maximum error over all processors
! ----------------------------------------------------------------------

        IF (lhook) CALL dr_hook('GCR_CALC_ABS_NORM',zhook_in,zhook_handle)
      i_start = 1
      i_end = row_length
      j_start = 1
      j_end = rows
      IF (model_domain  ==  mt_lam) THEN
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_end = rows-1
        IF(at_extremity(PEast)) i_end = row_length-1
        IF(at_extremity(PWest)) i_start = 2
      ELSE IF (model_domain  ==  mt_cyclic_lam) THEN
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_end = rows-1
      END IF

      Abs_norm= 0.0
      DO k = 1, model_levels
        DO j = j_start, j_end
          DO i = i_start, i_end
            IF (error(i,j,k)  >   abs_norm ) abs_norm = error(i,j,k)
          END DO
        END DO
      END DO

! Calculate max over all processors

      CALL gc_rmax(1, n_proc, info, abs_norm)

      IF (lhook) CALL dr_hook('GCR_CALC_ABS_NORM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_calc_abs_norm
