! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE Row_col_2norm

      SUBROUTINE row_col_2norm(                                         &
                               field, row_length, rows, levels,         &
                               level, halo_i, halo_j, off_u, off_v,     &
                               print_time, L_do_rims, rims_to_do )

      USE proc_info_mod
      USE timestep_mod

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

! Purpose:
!          Calculates row and column norms for 1 level of an array
!          Prints out every row norm and then every column norm
!          Output should use tkdiff for comparing runs.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      LOGICAL  L_do_rims ! LAMs include rims, Global do all polar points

      INTEGER                                                           &
        row_length,                                                     &
                         ! number of point on a row.
        rows,                                                           &
                         ! number of rows.
        levels,                                                         &
                         ! number of levels in field.
        level,                                                          &
                         ! level for norm calculation
        halo_i,                                                         &
        halo_j,                                                         &
        off_u,                                                          &
                         ! = 1 if field is at u-point
        off_v,                                                          &
                         ! = 1 if field is at v-point
        rims_to_do,                                                     &
        print_time       ! timestep_number for print output

      REAL, INTENT(IN) :: field(1-halo_i:row_length+halo_i,             &
                                1-halo_j:rows+halo_j, levels)

! Local Variables.

      INTEGER                                                           &
        i, j,                                                           &
        i_start,                                                        &
        i_stop,                                                         &
        j_start,                                                        &
        j_stop,                                                         &
        info,                                                           &
        gi,                                                             &
        gj,                                                             &
        istat

! Local arrays for parallel code

      REAL                                                              &
        row_wise( row_length, rows ),                                   &
        col_wise( rows, row_length ),                                   &
        two_norm_cols( row_length ),                                    &
        two_norm_rows( rows )

      REAL                                                              &
        global_2norm_rows(global_rows),                                 &
        global_2norm_cols(global_row_length)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1.   Calculate Norm.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ROW_COL_2NORM',zhook_in,zhook_handle)

      IF ( print_time == timestep_number ) THEN

        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows

        IF ( L_do_rims ) THEN

          IF ( model_domain == mt_global ) THEN
            IF(at_extremity(PNorth)) j_stop = rows - off_v
          END IF ! model_domain == mt_global

          IF( model_domain == mt_lam ) THEN
            IF(at_extremity(PNorth)) j_stop = rows - off_v
            IF(at_extremity(PEast)) i_stop = row_length - off_u
          END IF ! model_domain == mt_lam

        ELSE  ! do not do rims

          IF ( model_domain == mt_global ) THEN
            IF(at_extremity(PSouth)) j_start = 2
            IF(at_extremity(PNorth)) j_stop = rows - 1
          END IF ! model_domain == mt_global

          IF( model_domain == mt_lam ) THEN
            IF(at_extremity(PSouth)) j_start = 1 + rims_to_do
            IF(at_extremity(PNorth)) j_stop = rows - rims_to_do
            IF(at_extremity(PEast)) i_stop = row_length - rims_to_do
            IF(at_extremity(PWest)) i_start = 1 + rims_to_do
          ELSE IF (model_domain == mt_cyclic_LAM) THEN
            IF(at_extremity(PSouth)) j_start = 1 + rims_to_do
            IF(at_extremity(PNorth)) j_stop = rows - rims_to_do
          END IF ! model_domain == mt_lam

        END IF ! L_do_rims

        row_wise = 0.0
        col_wise = 0.0
        two_norm_rows = 0.0
        two_norm_cols = 0.0
        global_2norm_rows = 0.0
        global_2norm_cols = 0.0

        IF ( model_domain == mt_Global .and. .not. L_do_rims) THEN
! Global model only calculate norm for one of the polar points
!  if L_do_rims is .false.

          IF ( at_extremity(PSouth) .and. (datastart(1) == 1))          &
            row_wise(1,1) = row_wise(1,1) +                             &
                                field(1,1,level) * field(1,1,level)
          IF ( at_extremity(PNorth) .and. (datastart(1) == 1))          &
            row_wise(1,rows) = row_wise(1,rows) +                       &
                               field(1,rows,level) * field(1,rows,level)
        END IF  !  model_domain == mt_Global .and. .not. L_do_rims

        DO j = j_start , j_stop
          DO i = i_start, i_stop
            row_wise(i,j) = row_wise(i,j) +                             &
                                   field(i,j,level) * field(i,j,level)
          END DO
        END DO

        DO j = j_start , j_stop
          DO i = i_start, i_stop
            col_wise(j,i) = row_wise(i,j)
          END DO
        END DO

        CALL gcg_rvecsumr( rows , rows , 1, row_length, col_wise,       &
                           proc_col_group, istat, two_norm_cols )
        CALL gcg_rvecsumr( row_length , row_length , 1, rows, row_wise, &
                           proc_row_group, istat, two_norm_rows )

        DO j = 1, rows
          gj = datastart(2) + j - 1
          global_2norm_rows(gj) = two_norm_rows(j) / global_row_length
        END DO !  j = 1, rows
        DO i = 1, row_length
          gi = datastart(1) + i - 1
          global_2norm_cols(gi) = two_norm_cols(i) / global_rows
        END DO !  i = 1, row_length

        CALL gcg_rsumr(global_rows, proc_col_group, info,               &
                                                    global_2norm_rows)
        CALL gcg_rsumr(global_row_length, proc_row_group, info,         &
                                                    global_2norm_cols)

! Normalize by dividing by number of points

        DO j = 1, global_rows
          global_2norm_rows(j) = SQRT( global_2norm_rows(j) /           &
                                       global_row_length)
        END DO !  j = 1, global_rows
        DO i = 1, global_row_length
          global_2norm_cols(i) = SQRT(global_2norm_cols(i) /global_rows)
        END DO !  i = 1, global_row_length

        IF (me == 0) THEN
          WRITE(6,'('' rows '')')
          DO j = 1, global_rows, 4
            WRITE(6,'(''row '',I4,4E24.16)') j,                         &
                     ( global_2norm_rows(i), i = j, j+3 )
          END DO
          WRITE(6,'('' columns '')')
          DO j = 1, global_row_length, 4
            WRITE(6,'(''point '',I4,4E24.16)') j,                       &
                     ( global_2norm_cols(i), i = j, j+3 )
          END DO
        END IF !  me == 0

      END IF ! print_time == timestep_number

      IF (lhook) CALL dr_hook('ROW_COL_2NORM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE row_col_2norm
