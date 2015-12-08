! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Performs distribution of polar data from the semi-lagrangian
! interpolation routines

! Subroutine Interface:
SUBROUTINE interpolation_distribute_poles(                                 & 
                         sp_send, sp_levels, np_send, np_levels,           &
                         number_of_inputs, n_procx, ibase, proc_row_group, &
                         ime, dim_i_out, dim_j_out, dim_k_out,             &
                         at_extremity, l_high, l_mono, l_conserv,          &
                         data_out_mono, data_out_high)
                   

USE mpl, ONLY : &
    mpl_real

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParParams
IMPLICIT NONE

! Description:

!   Performs distribution of polar data from the semi-lagrangian
!   interpolation routines

! Method:
!   Uses a mpl_allgatherv to send data to all processors in row group

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

! Subroutine arguments

INTEGER, INTENT(IN) :: dim_i_out            ! x size of output arrays
INTEGER, INTENT(IN) :: dim_j_out            ! y size of output arrays
INTEGER, INTENT(IN) :: dim_k_out            ! z size of output arrays
INTEGER, INTENT(IN) :: ibase                ! base processor on row
INTEGER, INTENT(IN) :: ime                  ! My processor id
INTEGER, INTENT(IN) :: n_procx              ! number of procs in row group
INTEGER, INTENT(IN) :: number_of_inputs     ! 1, 2 or 3 inputs being interp.
INTEGER, INTENT(IN) :: proc_row_group       ! GCOM group for processor row

INTEGER, INTENT(IN) :: np_send  (0:n_procx-1)            ! Number of points 
                                                         ! sent at N pole
INTEGER, INTENT(IN) :: np_levels(0:n_procx-1, dim_k_out) ! levels sent at N Pole
INTEGER, INTENT(IN) :: sp_send  (0:n_procx-1)            ! Number of points 
                                                         ! sent at S pole
INTEGER, INTENT(IN) :: sp_levels(0:n_procx-1, dim_k_out) ! levels sent at S pole

LOGICAL, INTENT(IN) :: at_extremity(4)      ! At N,S,E,W extremity
LOGICAL, INTENT(IN) :: l_high               ! high order interpolation  = true
LOGICAL, INTENT(IN) :: l_mono               ! monotone interpolation = true
LOGICAL, INTENT(IN) :: l_conserv            ! monotone and conservative = true

REAL, INTENT(INOUT) :: data_out_mono(dim_i_out, dim_j_out, dim_k_out, &
                                     number_of_inputs)    ! data for polar
                                                          ! distribution
REAL, INTENT(INOUT) :: data_out_high(dim_i_out, dim_j_out, dim_k_out, &
                                     number_of_inputs)    ! data for polar
                                                          ! distribution

! Local variables

INTEGER :: info       ! GCOM errors
INTEGER :: i          ! Looper
INTEGER :: j          ! Looper
INTEGER :: k          ! Looper
INTEGER :: kk         ! Looper
INTEGER :: len        ! Message length
INTEGER :: n          ! Looper

INTEGER :: count (0:n_procx-1)   ! sizes to send
INTEGER :: offset(0:n_procx-1)   ! buffer offsets

REAL :: bcast_data(4*number_of_inputs*dim_k_out*n_procx) ! copy of data
                                                         ! for broadcast


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook('INTERPOLATION_DISTRIBUTE_POLES',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Section 1. South Pole Communications
!-----------------------------------------------------------------------
IF (at_extremity(psouth)) THEN

!-----------------------------------------------------------------------
! Section 1.1 South Pole Communications
! High order and monotone and conservative case
!-----------------------------------------------------------------------
  IF (l_high .AND. l_mono .AND. l_conserv) THEN

    offset(0) = 0
    count (0) = 2 * number_of_inputs * sp_send(0)

    DO j = 1, n_procx-1
      offset(j) = 2 * number_of_inputs * SUM(sp_send(0:j-1))
      count (j) = 2 * number_of_inputs * sp_send(j)
    END DO

    DO n = 1, number_of_inputs
      DO kk = 1, sp_send(ime)
        k = sp_levels(ime,kk)
        bcast_data(offset(ime) + 2*(n-1) * sp_send(ime) + 2*(kk-1) + 1) = &
                                       data_out_mono(1,1,k,n)
        bcast_data(offset(ime) + 2*(n-1) * sp_send(ime) + 2*(kk-1) + 2) = &
                                       data_out_high(1,1,k,n)
      END DO
    END DO

    CALL mpl_allgatherv(bcast_data(offset(ime)+1), count(ime),     mpl_real, &
                        bcast_data,                count,          offset,   &
                        mpl_real,                  proc_row_group, info )


    DO j = 0, n_procx-1
      IF (sp_send(j) > 0) THEN
        DO n = 1, number_of_inputs
          DO kk = 1, sp_send(j)
            k = sp_levels(j,kk)
            DO i = 1, dim_i_out
              data_out_mono(i,1,k,n) =                                &
                 bcast_data(offset(j) + 2*(n-1)*sp_send(j) + 2*(kk-1) + 1)
              data_out_high(i,1,k,n) =                                &
                 bcast_data(offset(j) + 2*(n-1)*sp_send(j) + 2*(kk-1) + 2)
            END DO
          END DO
        END DO
      END IF
    END DO

!-----------------------------------------------------------------------
! Section 1.2 South Pole Communications
! High order case
!-----------------------------------------------------------------------
  ELSE IF (l_high) THEN

    offset(0) = 0
    count (0) = number_of_inputs * sp_send(0)

    DO j = 1, n_procx-1
      offset(j) = number_of_inputs * SUM(sp_send(0:j-1))
      count (j) = number_of_inputs * sp_send(j)
    END DO

    DO n = 1, number_of_inputs
      DO kk = 1, sp_send(ime)
        k = sp_levels(ime,kk)
        bcast_data(offset(ime) + (n-1) * sp_send(ime) + kk) =        &
                                       data_out_high(1,1,k,n)
      END DO
    END DO

    CALL mpl_allgatherv(bcast_data(offset(ime)+1), count(ime),     mpl_real, &
                        bcast_data,                count,          offset,   &
                        mpl_real,                  proc_row_group, info )

    DO j = 0, n_procx-1
      IF (sp_send(j) > 0) THEN
        DO n = 1, number_of_inputs
          DO kk = 1, sp_send(j)
            k = sp_levels(j,kk)
            DO i = 1, dim_i_out
              data_out_high(i,1,k,n) =                                &
                 bcast_data(offset(j) + (n-1)*sp_send(j) + kk)
            END DO
          END DO
        END DO
      END IF
    END DO

!-----------------------------------------------------------------------
! Section 1.1 South Pole Communications
! Remainder of cases
!-----------------------------------------------------------------------
  ELSE

    offset(0) = 0
    count (0) = number_of_inputs * sp_send(0)

    DO j = 1, n_procx-1
      offset(j) = number_of_inputs * SUM(sp_send(0:j-1))
      count (j) = number_of_inputs * sp_send(j)
    END DO

    DO n = 1, number_of_inputs
      DO kk = 1, sp_send(ime)
        k = sp_levels(ime,kk)
        bcast_data(offset(ime) + (n-1) * sp_send(ime) + kk) =        &
                                       data_out_mono(1,1,k,n)
      END DO
    END DO

    CALL mpl_allgatherv(bcast_data(offset(ime)+1), count(ime),     mpl_real, &
                        bcast_data,                count,          offset,   &
                        mpl_real,                  proc_row_group, info )

    DO j = 0, n_procx-1
      IF (sp_send(j) > 0) THEN
        DO n = 1, number_of_inputs
          DO kk = 1, sp_send(j)
            k = sp_levels(j,kk)
            DO i = 1, dim_i_out
              data_out_mono(i,1,k,n) =                                &
                 bcast_data(offset(j) + (n-1)*sp_send(j) + kk)
            END DO
          END DO
        END DO
      END IF
    END DO

  END IF   ! l_high/mono/conserv

END IF     ! psouth

!-----------------------------------------------------------------------
! Section 2. North Pole Communications
!-----------------------------------------------------------------------
IF (at_extremity(pnorth)) THEN

!-----------------------------------------------------------------------
! Section 2.1 North Pole Communications
! High order and monotone and conservative case
!-----------------------------------------------------------------------
  IF (l_high .AND. l_mono .AND. l_conserv) THEN

    offset(0) = 0
    count (0) = 2 * number_of_inputs * np_send(0)

    DO j = 1, n_procx-1
      offset(j) = 2 * number_of_inputs * SUM(np_send(0:j-1))
      count (j) = 2 * number_of_inputs * np_send(j)
    END DO

    DO n = 1, number_of_inputs
      DO kk = 1, np_send(ime)
        k = np_levels(ime,kk)
        bcast_data(offset(ime) + 2*(n-1) * np_send(ime) + 2*(kk-1) + 1) = &
                                       data_out_mono(1,dim_j_out,k,n)
        bcast_data(offset(ime) + 2*(n-1) * np_send(ime) + 2*(kk-1) + 2) = &
                                       data_out_high(1,dim_j_out,k,n)
      END DO
    END DO

    CALL mpl_allgatherv(bcast_data(offset(ime)+1), count(ime),     mpl_real, &
                        bcast_data,                count,          offset,   &
                        mpl_real,                  proc_row_group, info )


    DO j = 0, n_procx-1
      IF (np_send(j) > 0) THEN
        DO n = 1, number_of_inputs
          DO kk = 1, np_send(j)
            k = np_levels(j,kk)
            DO i = 1, dim_i_out
              data_out_mono(i,dim_j_out,k,n) =                             &
                 bcast_data(offset(j) + 2*(n-1)*np_send(j) + 2*(kk-1) + 1)
              data_out_high(i,dim_j_out,k,n) =                             &
                 bcast_data(offset(j) + 2*(n-1)*np_send(j) + 2*(kk-1) + 2)
            END DO
          END DO
        END DO
      END IF
    END DO

!-----------------------------------------------------------------------
! Section 2.2 North Pole Communications
! High order case
!-----------------------------------------------------------------------
  ELSE IF (l_high) THEN

    offset(0) = 0
    count (0) = number_of_inputs * np_send(0)

    DO j = 1, n_procx-1
      offset(j) = number_of_inputs * SUM(np_send(0:j-1))
      count (j) = number_of_inputs * np_send(j)
    END DO

    DO n = 1, number_of_inputs
      DO kk = 1, np_send(ime)
        k = np_levels(ime,kk)
        bcast_data(offset(ime) + (n-1) * np_send(ime) + kk) =        &
                                       data_out_high(1,dim_j_out,k,n)
      END DO
    END DO

    CALL mpl_allgatherv(bcast_data(offset(ime)+1), count(ime),     mpl_real, &
                        bcast_data,                count,          offset,   &
                        mpl_real,                  proc_row_group, info )

    DO j = 0, n_procx-1
      IF (np_send(j) > 0) THEN
        DO n = 1, number_of_inputs
          DO kk = 1, np_send(j)
            k = np_levels(j,kk)
            DO i = 1, dim_i_out
              data_out_high(i,dim_j_out,k,n) =                         &
                 bcast_data(offset(j) + (n-1)*np_send(j) + kk)
            END DO
          END DO
        END DO
      END IF
    END DO

!-----------------------------------------------------------------------
! Section 2.1 North Pole Communications
! Remainder of cases
!-----------------------------------------------------------------------
  ELSE
    offset(0) = 0
    count (0) = number_of_inputs * np_send(0)

    DO j = 1, n_procx-1
      offset(j) = number_of_inputs * SUM(np_send(0:j-1))
      count (j) = number_of_inputs * np_send(j)
    END DO

    DO n = 1, number_of_inputs
      DO kk = 1, np_send(ime)
        k = np_levels(ime,kk)
        bcast_data(offset(ime) + (n-1) * np_send(ime) + kk) =        &
                                       data_out_mono(1,dim_j_out,k,n)
      END DO
    END DO

    CALL mpl_allgatherv(bcast_data(offset(ime)+1), count(ime),     mpl_real, &
                        bcast_data,                count,          offset,   &
                        mpl_real,                  proc_row_group, info )

    DO j = 0, n_procx-1
      IF (np_send(j) > 0) THEN
        DO n = 1, number_of_inputs
          DO kk = 1, np_send(j)
            k = np_levels(j,kk)
            DO i = 1, dim_i_out
              data_out_mono(i,dim_j_out,k,n) =                         &
                 bcast_data(offset(j) + (n-1)*np_send(j) + kk)
            END DO
          END DO
        END DO
      END IF
    END DO

  END IF !   l_high/mono/conserv

END IF   ! pnorth

IF (lhook) CALL dr_hook('INTERPOLATION_DISTRIBUTE_POLES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE interpolation_distribute_poles
