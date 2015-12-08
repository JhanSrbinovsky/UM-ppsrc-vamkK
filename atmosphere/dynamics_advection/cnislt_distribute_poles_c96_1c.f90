! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Performs distribution of polar data from the calc_non_int_sl_theta

! Subroutine Interface:
SUBROUTINE cnislt_distribute_poles(                              &
               sp_send, sp_levels, np_send, np_levels,           &
               n_procx, ibase, proc_row_group,                   &
               rows, row_length, model_levels,                   & 
               off_x1, off_y1, off_x2, off_y2,                   &
               ime,  at_extremity, w_logic,                      &
               field_1, field_2)


USE mpl, ONLY : &
    mpl_real

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParParams
IMPLICIT NONE

! Description:

!   Performs distribution of polar data from calc_non_int_sl_theta

! Method:
!   Uses a loop around  processors in polar row doing broadcasts

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

! Subroutine arguments

INTEGER, INTENT(IN) :: ibase                ! base processor on row
INTEGER, INTENT(IN) :: ime                  ! My processor id
INTEGER, INTENT(IN) :: model_levels         ! Number of levels
INTEGER, INTENT(IN) :: n_procx              ! number of procs in row group
INTEGER, INTENT(IN) :: off_x1               ! Field 1 x halo
INTEGER, INTENT(IN) :: off_y1               ! Field 1 y halo
INTEGER, INTENT(IN) :: off_x2               ! Field 2 x halo
INTEGER, INTENT(IN) :: off_y2               ! Field 2 y halo
INTEGER, INTENT(IN) :: proc_row_group       ! GCOM group for processor row
INTEGER, INTENT(IN) :: row_length           ! Length of row in fields
INTEGER, INTENT(IN) :: rows                 ! Number of rows in fields

INTEGER, INTENT(IN) :: np_send  (0:n_procx-1)            ! Number of points
                                                         ! sent at N pole
INTEGER, INTENT(IN) :: np_levels(0:n_procx-1, model_levels) ! levels sent at 
                                                            ! N Pole
INTEGER, INTENT(IN) :: sp_send  (0:n_procx-1)            ! Number of points
                                                         ! sent at S pole
INTEGER, INTENT(IN) :: sp_levels(0:n_procx-1, model_levels) ! levels sent at
                                                            ! S pole

LOGICAL, INTENT(IN) :: w_logic              ! use the logic for w fields
LOGICAL, INTENT(IN) :: at_extremity(4)      ! At N,S,E,W extremity

REAL, INTENT(INOUT) :: field_1(1-off_x1:row_length+off_x1, &
                               1-off_y1:rows+off_y1,       &
                               model_levels) ! data for polar distribution
REAL, INTENT(INOUT) :: field_2(1-off_x2:row_length+off_x2, &
                               1-off_y2:rows+off_y2,       &
                               model_levels) ! data for polar distribution

! Local variables

INTEGER :: info       ! GCOM errors
INTEGER :: i          ! Looper
INTEGER :: j          ! Looper
INTEGER :: k          ! Looper
INTEGER :: kk         ! Looper

INTEGER :: offset(0:n_procx-1)          ! Offsets for gatherv
INTEGER :: count (0:n_procx-1)          ! Counts for gatherv

REAL :: ctmp1(2*model_levels*n_procx)   ! Copy of data for broadcast
REAL :: ctmp2(2,model_levels)   ! Copy of data for broadcast


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook('CNISLT_DISTRIBUTE_POLES',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Section 1. South Pole Communications
!-----------------------------------------------------------------------
IF (at_extremity(psouth)) THEN

! Setup offsets and counts
  offset(0) = 0
  count (0) = 2 * sp_send(0)

  DO j = 1, n_procx-1
    count (j) = 2 * sp_send(j)
    offset(j) = 2 * SUM(sp_send(0:j-1))
  END DO

! Put data into buffer for comms
  DO kk = 1, sp_send(ime)
    k = sp_levels(ime,kk)
    ctmp1(offset(ime) + 2*kk -1) = field_1(1,1,k)
    ctmp1(offset(ime) + 2*kk   ) = field_2(1,1,k)
  END DO

  CALL mpl_allgatherv(ctmp1(offset(ime)+1), count(ime),     mpl_real,   &
                      ctmp1,                count,          offset,     &
                      mpl_real,             proc_row_group, info)


! Unpack the data
  DO j = 0, n_procx-1
    IF (sp_send(j) > 0) THEN
      DO kk = 1, sp_send(j)
        k = sp_levels(j,kk)

! W/W_star need different logic from the rest of the fields
        IF (w_logic) THEN
          IF (ctmp1(offset(j) + 2*kk - 1)  ==                        & 
              ctmp1(offset(j) + 2*kk    ))    THEN
            DO i = 1, row_length
              field_2(i,1,k) = field_1(i,1,k)
            END DO
          ELSE
            DO i = 1, row_length
              field_2(i,1,k) = ctmp1(offset(j) + 2*kk)
            END DO
          END IF

        ELSE
          DO i = 1, row_length
            field_1(i,1,k) = ctmp1(offset(j) + 2*kk -1)
            field_2(i,1,k) = ctmp1(offset(j) + 2*kk   )
          END DO
        END IF

      END DO
    END IF    ! sp_send
  END DO

END IF

!-----------------------------------------------------------------------
! Section 2. North Pole Communications
!-----------------------------------------------------------------------
IF (at_extremity(pnorth)) THEN

! Setup offsets and counts
  offset(0) = 0
  count (0) = 2 * np_send(0)

  DO j = 1, n_procx-1
    count (j) = 2 * np_send(j)
    offset(j) = 2 * SUM(np_send(0:j-1))
  END DO

! Put data in buffer for comms
  DO kk = 1, np_send(ime)
    k = np_levels(ime,kk)
    ctmp1(offset(ime) + 2*kk -1) = field_1(1,rows,k)
    ctmp1(offset(ime) + 2*kk   ) = field_2(1,rows,k)
  END DO

  CALL mpl_allgatherv(ctmp1(offset(ime)+1), count(ime),     mpl_real,   &
                      ctmp1,                count,          offset,     &
                      mpl_real,             proc_row_group, info)


  DO j = 0, n_procx-1
    IF (np_send(j) > 0) THEN
      DO kk = 1, np_send(j)
        k = np_levels(j,kk)

! W/W_star need different logic from the rest of the fields
        IF (w_logic) THEN
          IF (ctmp1(offset(j) + 2*kk - 1)  ==                        & 
              ctmp1(offset(j) + 2*kk    ))    THEN
            DO i = 1, row_length
              field_2(i,rows,k) = field_1(i,rows,k)
            END DO
          ELSE
            DO i = 1, row_length
              field_2(i,rows,k) = ctmp1(offset(j) + 2*kk)
            END DO
          END IF

        ELSE
          DO i = 1, row_length
            field_1(i,rows,k) = ctmp1(offset(j) + 2*kk -1)
            field_2(i,rows,k) = ctmp1(offset(j) + 2*kk   )
          END DO
        END IF

      END DO
    END IF    ! np_send
  END DO

END IF

IF (lhook) CALL dr_hook('CNISLT_DISTRIBUTE_POLES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE cnislt_distribute_poles
