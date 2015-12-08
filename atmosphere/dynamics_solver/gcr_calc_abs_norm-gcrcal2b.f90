! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine GCR_calc_abs_norm_2B
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

      SUBROUTINE GCR_calc_abs_norm_2B(                                  &
     &                            Error, off_x, off_y,                  &
     &                            i_start, i_stop, j_start, j_stop,     &
     &                            n_proc, row_length, rows,             &
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
     &, i_start, i_stop                                                 &
                                   ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop            ! loop bounds set in PE_Helmholtz

      REAL                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Arguments with Intent OUT. ie: variables Output only
      REAL                                                              &
     &  Abs_norm

! Local variables

      INTEGER i, j, k, info

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External Routines:

      EXTERNAL                                                          &
     &  gc_rmax

! ----------------------------------------------------------------------
! Section 1.   Calculate maximum error over all processors
! ----------------------------------------------------------------------


      IF (lhook) CALL dr_hook('GCR_CALC_ABS_NORM_2B',zhook_in,zhook_handle)
      Abs_norm= 0.0
      DO k = 1, model_levels
        DO j = j_start, j_stop
          DO i = i_start, i_stop
            IF (error(i,j,k)  >   abs_norm ) abs_norm = error(i,j,k)
          END DO
        END DO
      END DO

! Calculate max over all processors

      CALL gc_rmax(1, n_proc, info, abs_norm)

      IF (lhook) CALL dr_hook('GCR_CALC_ABS_NORM_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_calc_abs_norm_2B
