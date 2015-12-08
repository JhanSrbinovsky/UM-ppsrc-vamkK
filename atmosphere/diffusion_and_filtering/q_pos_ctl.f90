! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine  Q_Pos_Ctl

SUBROUTINE q_pos_ctl                                              &
          (q, row_length, rows, wet_model_levels, variables,      &
           bl_levels, global_row_length, global_rows,             &
           me, n_proc, off_x, off_y,                              &
           model_domain,                                          &
           halo_type, q_pos_method, qlimit,                       &
           l_qpos_diag_pr, qpos_diag_limit, ident)


! Method:
!         Multiple methods for resetting Q are now supported.
!         All methods (except original non-local) will invoke a polar
!         row reset for global models.
!         For reset, level, column and hybrid Q_POS methods simple
!         subroutine calls deal with the data in place. For the original
!         and local methods the data is first gathered level-by-level
!         onto processors where the algorithms are applied before the
!         results are returned to their original processors. These latter
!         two methods are much more costly and much less scalable as a result.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

USE q_pos_method_mod, ONLY : q_pos_original, q_pos_local, q_pos_reset,  &
                             q_pos_column,   q_pos_level, q_pos_hybrid

USE global_2d_sums_mod, ONLY : global_2d_sums

USE dynamics_input_mod, ONLY: l_endgame

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY : ereport
USE UM_ParVars
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) ::  row_length          ! number of points on a row
INTEGER, INTENT(IN) ::  rows                ! number of rows of data
INTEGER, INTENT(IN) ::  wet_model_levels    ! number of levels used
INTEGER, INTENT(IN) ::  variables           ! number of variables
INTEGER, INTENT(IN) ::  bl_levels           ! number of boundary layer levels
INTEGER, INTENT(IN) ::  global_row_length   ! num of points in row in full field
INTEGER, INTENT(IN) ::  global_rows         ! number of rows in full field
INTEGER, INTENT(IN) ::  me                  ! processor number
INTEGER, INTENT(IN) ::  n_proc              ! number of processors
INTEGER, INTENT(IN) ::  off_x               ! Size of halo in i direction.
INTEGER, INTENT(IN) ::  off_y               ! Size of halo in j direction.
INTEGER, INTENT(IN) ::  halo_type           ! halo type
INTEGER, INTENT(IN) ::  model_domain        ! model domain type
INTEGER, INTENT(IN) ::  q_pos_method        ! Method used

LOGICAL, INTENT(IN) :: l_qpos_diag_pr    ! Are diagnostic prints enabled?

REAL, INTENT(IN)    :: qlimit            ! Limit for Q
REAL, INTENT(IN)    :: qpos_diag_limit   ! Limit for diagnostic prints

CHARACTER (LEN=*), INTENT(IN) :: ident   ! Identify which q_pos call this is

! Arguments with Intent IN/OUT.

REAL, INTENT(INOUT) ::  q(1-off_x:row_length+off_x, 1-off_y:rows+off_y, &
                          wet_model_levels)


! Local Variables.
INTEGER :: i, j, k                     ! Loopers
INTEGER :: max_levs_per_cpu            ! max number of levs per cpu
INTEGER :: icode                       ! return code from comms.
INTEGER :: proc                        ! processor being sent to
INTEGER :: actual_levels
LOGICAL :: l_q_pos_local               ! for method 1/2

REAL :: recip_g_row_len                ! 1/global_row_length
REAL :: polar_value                    ! new constant value at pole

INTEGER :: map(wet_model_levels)       ! processor number for level
INTEGER :: n_levs_on_proc(0:n_proc-1)  ! number of levels on each processor
INTEGER :: levels_to_fix(wet_model_levels) ! Which levels need fixing?

REAL :: local_polar(row_length, wet_model_levels) ! local copy of pole
REAL :: polar_sum(wet_model_levels)    ! Summation of polar values
REAL ::  q_global (global_row_length, global_rows,                      &
                   (wet_model_levels/n_proc)+1) ! Gathered Q


CHARACTER (LEN=80)           :: cmessage         ! error message
CHARACTER (LEN=*), PARAMETER :: routinename = 'Q_POS_CTL'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!========================================================================
IF (lhook) CALL dr_hook('Q_POS_CTL',zhook_in,zhook_handle)

!-------------------------------------------------------------------------
! Section 0. Initialise options
!-------------------------------------------------------------------------
IF (q_pos_method == q_pos_original) THEN
  l_q_pos_local = .FALSE.
ELSE IF (q_pos_method == q_pos_local) THEN
  l_q_pos_local = .TRUE.
END IF

!-------------------------------------------------------------------------
! Section 0.1. Perform diagnostics if required
!-------------------------------------------------------------------------
IF (l_qpos_diag_pr) THEN
! DEPENDS ON: q_pos_diag_pr
  CALL q_pos_diag_pr(q, row_length, rows, wet_model_levels,             &
                     off_x, off_y, qpos_diag_limit, ident)
END IF

!-------------------------------------------------------------------------
! Section 1. Reset Polar rows to constant for global model if using
!            local q_pos method or reset/column/level
!-------------------------------------------------------------------------
IF ( model_domain == mt_global    .AND.                                 &
     .NOT. l_endgame              .AND.                                 &
    (q_pos_method == q_pos_local  .OR.                                  &
     q_pos_method == q_pos_reset  .OR.                                  &
     q_pos_method == q_pos_column .OR.                                  &
     q_pos_method == q_pos_level  .OR.                                  &
     q_pos_method == q_pos_hybrid)       ) THEN

  IF (at_extremity(pnorth)) THEN

    DO k = 1, wet_model_levels
      DO i = 1, row_length
        local_polar(i, k) = q(i, rows, k)
      END DO
    END DO

    CALL global_2d_sums(local_polar, row_length, 1, 0, 0,               &
                        wet_model_levels, polar_sum,                    &
                        gc_proc_row_group)

    recip_g_row_len = 1./global_row_length
    DO k = 1, wet_model_levels
      polar_value   = polar_sum(k) * recip_g_row_len
      q(:, rows, k) = polar_value
    END DO

  END IF  ! PNORTH

  IF (at_extremity(psouth)) THEN

    DO k = 1, wet_model_levels
      DO i = 1, row_length
        local_polar(i, k) = q(i, 1, k)
      END DO
    END DO

    CALL global_2d_sums(local_polar, row_length, 1, 0, 0,               &
                        wet_model_levels, polar_sum,                    &
                        gc_proc_row_group)

    recip_g_row_len = 1./global_row_length
    DO k = 1, wet_model_levels
      polar_value   = polar_sum(k) * recip_g_row_len
      q(:, 1, k) = polar_value
    END DO

  END IF  ! PSOUTH
END IF ! GLOBAL and not endgame

IF (q_pos_method == q_pos_reset) THEN
! DEPENDS ON : q_pos_reset_sub
  CALL q_pos_reset_sub(q, row_length, rows, off_x, off_y,               &
                       wet_model_levels, qlimit )

ELSE IF (q_pos_method == q_pos_column) THEN
  actual_levels = wet_model_levels/variables
! DEPENDS ON : q_pos_column_sub
  DO i = 1, variables
    CALL q_pos_column_sub(q(:,:,(i-1)*actual_levels+1:i*actual_levels), &
                          row_length, rows, off_x, off_y,               &
                          actual_levels, qlimit )
  END DO

ELSE IF (q_pos_method == q_pos_level) THEN
! DEPENDS ON : q_pos_level_sub
  CALL q_pos_level_sub(q, row_length, rows, off_x, off_y,               &
                       wet_model_levels, qlimit )

ELSE IF (q_pos_method == q_pos_hybrid) THEN
! DEPENDS ON : q_pos_hybrid_sub
  CALL q_pos_hybrid_sub(q, row_length, rows, off_x, off_y,              &
                       wet_model_levels, variables, bl_levels, qlimit )

ELSE IF (q_pos_method == q_pos_original .OR.                            &
         q_pos_method == q_pos_local)  THEN
!-------------------------------------------------------------------------
! Section 2.  Check all levels to see which ones have values locally
!             that need correcting and setup maps
!-------------------------------------------------------------------------
  levels_to_fix(:) = 0

  DO k=1, wet_model_levels
    outer: DO j=1, rows
             DO i=1, row_length
               IF ( q(i,j,k) < qlimit ) THEN
                 levels_to_fix(k) = 1

                 ! We have found a point that needs fixing, so can jump
                 ! out of the loop.  This would inhibit vectorisation
                 ! so only do it for non-vector platforms
                 EXIT outer

               END IF
             END DO
    END DO outer
  END DO

  ! We need to get the maximum of this across all processors to find
  ! those levels that need fixing
  CALL gc_imax(wet_model_levels, n_proc, icode, levels_to_fix)

  ! Calculate the mapping of which processor each level will go to
  n_levs_on_proc(:) = 0
  map(:)            = -99
  proc              = 0

  DO k=1, wet_model_levels
    IF (levels_to_fix(k) > 0) THEN
      map(k)               = proc
      proc                 = proc + 1
      IF (proc == n_proc) THEN
       proc                = 0
      END IF
    END IF
  END DO

!-------------------------------------------------------------------------
! Section 3.  Distribute Q over the processors
!-------------------------------------------------------------------------

  DO k=1, wet_model_levels
    IF (levels_to_fix(k) > 0) THEN

      n_levs_on_proc(map(k)) = n_levs_on_proc(map(k)) + 1

! DEPENDS ON: gather_field
      CALL gather_field( q(1-off_x,1-off_y,k),                      &
                         q_global(1,1,n_levs_on_proc(map(k))),      &
                         row_length + 2*off_x, rows + 2*off_y,      &
                         global_row_length, global_rows,            &
                         fld_type_p, halo_type,                     &
                         map(k), gc_all_proc_group,                 &
                         icode, cmessage)

    END IF
  END DO

!-------------------------------------------------------------------------
! Section 4.  Do QPOS for all the levels on this processor
!-------------------------------------------------------------------------

! DEPENDS ON: q_pos
  CALL q_pos(                                                       &
             q_global, global_row_length, global_rows,              &
             n_levs_on_proc(me),                                    &
             l_q_pos_local, qlimit ,                                &
             model_domain)


!-------------------------------------------------------------------------
! Section 5.  Scatter back levels to original processors
!-------------------------------------------------------------------------
  n_levs_on_proc(:) = 0
  DO k=1, wet_model_levels
    IF (levels_to_fix(k) > 0) THEN

      n_levs_on_proc(map(k)) = n_levs_on_proc(map(k)) + 1
! DEPENDS ON: scatter_field
      CALL scatter_field( q(1-off_x,1-off_y,k),                     &
                          q_global(1,1,n_levs_on_proc(map(k))),     &
                          row_length + 2*off_x, rows + 2*off_y,     &
                          global_row_length, global_rows,           &
                          fld_type_p, halo_type,                    &
                          map(k), gc_all_proc_group,                &
                          icode, cmessage)

    END IF
  END DO

ELSE ! no q_pos_method I recognise
  icode = -10
  cmessage = 'Q_POS: unrecognised method'

  CALL ereport(routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook('Q_POS_CTL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE q_pos_ctl

