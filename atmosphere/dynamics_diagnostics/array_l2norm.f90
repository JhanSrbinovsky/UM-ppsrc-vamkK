! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine array_l2norm

      SUBROUTINE array_l2norm(                                          &
                              field, row_length, rows, levels,          &
                              start_level, end_level,                   &
                              halo_i, halo_j, off_u, off_v,             &
                              L_do_halos, L_do_rims,                    &
                              rim_ew, rim_ns, Two_Norm )

      USE proc_info_mod

! Purpose:
!          Calculates the two norm over a range of levels of the
!          input field.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE global_2d_sums_mod, ONLY: &
          global_2d_sums,           &
          global_sum_fast,          &
          global_sum_reprod_dd,     &
          global_sum_method

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      LOGICAL                                                           &
        L_do_rims,                                                      &
                     ! LAMs: include rims, Global do all polar points
        L_do_halos     ! include halos in calculation of norm
!                    NB this means that some points are counted twice

      INTEGER                                                           &
        row_length,                                                     &
                       ! number of point on a row.
        rows,                                                           &
                       ! number of rows.
        levels,                                                         &
                       ! number of levels in field.
        start_level,                                                    &
                       ! start_level for norm calculation
        end_level      ! end_level for norm calculation

      INTEGER                                                           &
        halo_i,                                                         &
        halo_j,                                                         &
        off_u,                                                          &
                       ! = 1 if field is at u-point
        off_v,                                                          &
                       ! = 1 if field is at v-point
        rim_ew,                                                         &
        rim_ns         ! rim specified for each direction since
                       ! solver coeff arrays span different points

      REAL, INTENT(IN) :: field(1-halo_i:row_length+halo_i,             &
                                1-halo_j:rows+halo_j, levels)

! Arguments with Intent OUT. ie: variables Output only

      REAL, INTENT(OUT) :: Two_Norm

! Local Variables.

      INTEGER                                                           &
        i, j, k,                                                        &
        i_start,                                                        &
        i_stop,                                                         &
        j_start,                                                        &
        j_stop,                                                         &
        ipoints,                                                        &
        istat

      INTEGER :: sum_method_old     ! temporary store for existing method
      REAL :: points     ! number of points used in norm calculation

! Local arrays for parallel code

      REAL :: two_norm_component( row_length+2*halo_i, rows+2*halo_j )
      REAL :: two_norm_sum( 1 )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1.   Calculate Norm.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ARRAY_L2NORM',zhook_in,zhook_handle)

      i_start = 1
      i_stop = row_length
      j_start = 1
      j_stop = rows

      IF ( L_do_rims ) THEN

        ipoints = global_row_length * global_rows

        IF( L_do_halos ) THEN
          i_start = 1 - halo_i
          i_stop = row_length + halo_i
          j_start = 1 - halo_j
          j_stop = rows + halo_j
          ipoints = (global_row_length + 2 * n_procx * halo_i ) *       &
                    (global_rows + 2 * n_procy * halo_j )
        END IF ! L_do_halos

        IF ( model_domain == mt_global ) THEN
          IF(at_extremity(PNorth)) j_stop = rows - off_v
          ipoints = global_row_length * (global_rows - off_v)
        END IF ! model_domain == mt_global

        IF( model_domain == mt_lam ) THEN
          IF(at_extremity(PNorth)) j_stop = rows - off_v
          IF(at_extremity(PEast)) i_stop = row_length - off_u
          ipoints = (global_row_length - off_u) * (global_rows - off_v)
        END IF ! model_domain == mt_lam

      ELSE  ! do not do halos or rims

        IF( L_do_halos ) THEN
          write(6,*)' You must do rims if halos are included in L2norms'
        END IF ! L_do_halos

        IF ( model_domain == mt_global ) THEN
          IF(at_extremity(PSouth)) j_start = 2
          IF(at_extremity(PNorth)) j_stop = rows - 1
          ipoints = global_row_length * (global_rows - 2) + 2
        END IF ! model_domain == mt_global

        IF( model_domain == mt_lam ) THEN
          IF(at_extremity(PSouth)) j_start = 1 + rim_ns
          IF(at_extremity(PNorth)) j_stop = rows - rim_ns - off_v
          IF(at_extremity(PEast)) i_stop = row_length - rim_ew - off_u
          IF(at_extremity(PWest)) i_start = 1 + rim_ew
          ipoints = (global_row_length - 2*rim_ew - off_u) *            &
                    (global_rows  - 2*rim_ns - off_v)

        ELSE IF (model_domain == mt_cyclic_LAM) THEN
          IF(at_extremity(PSouth)) j_start = 1 + rim_ns
          IF(at_extremity(PNorth)) j_stop = rows - rim_ns
          ipoints = global_row_length * (global_rows - 2*rim_ns)

        ELSE IF (model_domain == mt_bi_cyclic_lam) THEN
          
          ipoints = global_row_length * global_rows

        END IF ! model_domain == mt_lam

      END IF ! L_do_rims

      points = REAL(ipoints)

      two_norm_component(:,:) = 0.0

      IF ( model_domain == mt_Global .AND. .NOT. L_do_rims) THEN
! Global model only calculate norm for one of the polar points
!  IF L_do_rims is .false.

        IF(at_extremity(PSouth) .AND. (l_datastart(1) == 1)) THEN
          DO k = start_level, end_level
            two_norm_component(1,1)= two_norm_component(1,1) +          &
                                     field(1,1,k) * field(1,1,k)
          END DO
        END IF
        IF(at_extremity(PNorth) .AND. (l_datastart(1) == 1)) THEN
          DO k = start_level, end_level
            two_norm_component(1,rows)= two_norm_component(1,rows) +    &
                                   field(1,rows,k) * field(1,rows,k)
          END DO
        END IF
      END IF  !  model_domain == mt_Global .AND. .NOT. L_do_rims

      DO k = start_level, end_level
        DO j = j_start , j_stop
          DO i = i_start, i_stop
            two_norm_component( i + halo_i, j + halo_j ) =              &
                two_norm_component( i + halo_i, j + halo_j ) +          &
                                      field(i,j,k) * field(i,j,k)
          END DO
        END DO
      END DO  ! k = start_level, end_level

      ! Whatever options chosen in the UMUI we want to force this sum to
      ! be reproducible as this is a debugging routine! Thus store old
      ! method and restore it after summation.
      sum_method_old = global_sum_method
      IF (global_sum_method == global_sum_fast) THEN
        global_sum_method = global_sum_reprod_dd
      END IF

      CALL global_2d_sums(two_norm_component, row_length + 2*halo_i,   &
                          rows + 2*halo_j, 0, 0, 1, two_norm_sum)

      global_sum_method = sum_method_old

! Normalize by dividing by number of points

      Two_Norm = SQRT (Two_Norm_Sum(1) / points)

      IF (lhook) CALL dr_hook('ARRAY_L2NORM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE array_l2norm
