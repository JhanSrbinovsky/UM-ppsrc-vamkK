! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module GLOBAL_2D_SUMS_MOD
MODULE global_2d_sums_mod


! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Purpose:
!    This module contains a number of methods for performing global
!    2D (and by limitation of dimension 1D) sums across a group of 
!    processors.
!
! Code Description:
!   Language: FORTRAN 90

USE UM_ParVars
IMPLICIT NONE

! Parameters for the different possible methods,
INTEGER, PARAMETER :: global_sum_reprod_orig = 1  ! Original method
INTEGER, PARAMETER :: global_sum_reprod_dd   = 2  ! "double double" method
                                                  ! in GCOM
INTEGER, PARAMETER :: global_sum_fast        = 3  ! Fast non-reprod method

INTEGER :: global_sum_method = global_sum_reprod_orig  ! Set a default

CONTAINS
  SUBROUTINE global_2d_sums(field, row_length, rows, off_x, off_y, &
                            levels, global_sum, gid)

! Purpose: this routine provides the global sums according to some 
!          the options chosen. Group ID is an optional argument but
!          is ignored for the original method.


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

  ! Arguments 
  INTEGER, INTENT(IN)  :: row_length
  INTEGER, INTENT(IN)  :: rows
  INTEGER, INTENT(IN)  :: off_x
  INTEGER, INTENT(IN)  :: off_y
  INTEGER, INTENT(IN)  :: levels
  INTEGER, INTENT(IN), OPTIONAL  :: gid   ! Group ID 

  REAL,    INTENT(IN)  :: field(1-off_x : row_length+off_x,  &
                                1-off_y : rows+off_y,        &
                                levels)

  REAL,    INTENT(OUT) :: global_sum(levels)

  ! Local variables
  INTEGER                       :: i,j,k
  INTEGER                       :: istat
  INTEGER                       :: my_gid
  REAL                          :: local_sum(MAX(levels,rows))
  CHARACTER (LEN=256)           :: cmessage
  CHARACTER (LEN=*), PARAMETER :: routinename='GLOBAL_2D_SUMS'

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
  IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)

  ! Deal with optional group
  IF (PRESENT(gid)) THEN
    my_gid = gid
  ELSE
    ! Default position is to use the current global communicator
    CALL gc_get_communicator(my_gid, istat)
    IF (istat /= 0) THEN
      WRITE(cmessage,*) 'Problem in gc_get_communicator'
      CALL ereport(routinename, 50, cmessage  )
    END IF

  END IF


  SELECT CASE (global_sum_method)
    CASE (global_sum_reprod_orig)

      ! Loop over levels since we could have an offset for rows so summation is
      ! not in memory order which GCOM expects since the offset is between
      ! levels.
      DO k = 1, levels

        ! Sum first along rows
        CALL gcg_rvecsumr(row_length+2*off_x, row_length, 1+off_x,       &
                          rows, field(:,1:,k), gc_proc_row_group, &
                          istat, local_sum)
        IF (istat /= 0) THEN
          WRITE(cmessage,*) 'Problem in gcg_rvecsumr'
          CALL ereport(routinename, 100, cmessage  )
        END IF
  
        ! If we have more than one row, sum up the rows
        IF (rows > 1) THEN
          CALL gcg_rvecsumr(rows, rows, 1, 1, local_sum,            &
                            gc_proc_col_group, istat, global_sum(k))
          IF (istat /= 0) THEN
            WRITE(cmessage,*) 'Problem in 2nd gcg_rvecsumr'
            CALL ereport(routinename, 200, cmessage  )
          END IF
        ELSE
          global_sum(k) = local_sum(1)
        END IF
      END DO

    CASE (global_sum_reprod_dd)
      ! Use the GCOM routine to do this
      CALL gcg_r2darrsum(field, row_length, rows, off_x, off_y, levels, &
                         my_gid, global_sum, istat)
      IF (istat /= 0) THEN
        WRITE(cmessage,*) 'Problem in gcg_r2darrsum'
        CALL ereport(routinename, 300, cmessage  )
      END IF

    CASE (global_sum_fast)
      ! Create local partial sums
      DO k = 1, levels
        local_sum(k) = 0.0

        DO j = 1, rows
          DO i = 1, row_length
            local_sum(k) = local_sum(k) + field(i,j,k)
          END DO
        END DO
      END DO

      ! Sum up partial sums
      CALL GCG_RSUM(levels, my_gid, istat, local_sum)
      IF (istat /= 0) THEN
        WRITE(cmessage,*) 'Problem in gcg_rsum'
        CALL ereport(routinename, 400, cmessage  )
      END IF

      global_sum(1:levels) = local_sum(1:levels)

    CASE DEFAULT
      WRITE(cmessage,*) 'Unrecognised global 2D sum method ', &
                         global_sum_method
      CALL ereport(routinename, 500, cmessage  )

  END SELECT

  IF (lhook) CALL dr_hook(routinename,zhook_out,zhook_handle)
  RETURN
  END SUBROUTINE global_2d_sums

END MODULE global_2d_sums_mod
