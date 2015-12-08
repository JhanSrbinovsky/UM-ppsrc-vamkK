! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_Vector_Wind

      Subroutine Polar_Vector_Wind_n(                                   &
                                     uv_comp,                           &
                                     sin_longitude, cos_longitude,      &
                                     row_length, in_rows, model_levels, &
                                     mag_vector_np, dir_vector_np,      &
                                     mag_vector_sp, dir_vector_sp,      &
                                     halo_i, halo_j,                    &
                                     global_row_length,                 &
                                     proc_row_group, at_extremity)

      USE global_2d_sums_mod, ONLY: global_2d_sums


! Purpose:
!          Calculates vector wind at each pole and returns magnitude
!          and direction.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      INTEGER                                                           &
        row_length,                                                     &
                     ! Dimension of input wind  in i direction.
        in_rows,                                                        &
                     ! Dimension of input wind in j direction
        model_levels,                                                   &
                     ! Dimension of  input wind in k direction.
        halo_i,                                                         &
                     ! Size of halo in i direction.
        halo_j,                                                         &
                     ! Size of halo in j direction.
        global_row_length,                                              &
                     ! global number of points on a row
        proc_row_group ! Group id for processors on the same row

      REAL :: uv_comp (1-halo_i:row_length+halo_i,                      &
                       1-halo_j:in_rows+halo_j, model_levels)
                                                 ! input wind component
      REAL                                                              &
        cos_longitude(row_length,in_rows),                              &
                                          ! cosine of longitude
        sin_longitude(row_length,in_rows) ! sine of longitude

      LOGICAL ::  at_extremity(4)  ! Indicates IF this processor is at 
                                   ! north ,south, east or west of domain

! Arguments with Intent OUT. ie: Output variables.

      REAL                                                              &
        mag_vector_np(model_levels),                                    &
                                    ! mag. of the vector wind at N pole
        dir_vector_np(model_levels),                                    &
                                  ! dirn. of the vector wind at N pole
        mag_vector_sp(model_levels),                                    &
                                  ! mag. of the vector wind at S pole
        dir_vector_sp(model_levels) ! dirn. of the vector wind at S pole

! Local Variables.

      INTEGER                                                           &
        i, k,                                                           &
               ! Loop indices
        info   ! Status variable

      REAL ::  a_np, b_np, a_sp, b_sp

      REAL ::  wrk(2*model_levels), rwrk(row_length,2*model_levels)

! Parameters
      INTEGER                                                           &
        PNorth,                                                         &
                      ! North processor address in the neighbor array
        PEast,                                                          &
                      ! East processor address in the neighbor array
        PSouth,                                                         &
                      ! South processor address in the neighbor array
        PWest,                                                          &
                      ! West processor address in the neighbor array
        NoDomain      ! Value in neighbor array IF the DOmain has
                       !  no neighbor in this direction. Otherwise
                       !  the value will be the tid of the neighbor
      PARAMETER (                                                       &
        PNorth   = 1,                                                   &
        PEast    = 2,                                                   &
        PSouth   = 3,                                                   &
        PWest    = 4,                                                   &
        NoDomain = -1)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Calculate magnitude and direction of polar vector wind
!              from v component of wind on row around pole.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('POLAR_VECTOR_WIND_N',zhook_in,zhook_handle)

      IF (at_extremity(PNorth)) THEN

         DO k = 1, model_levels
            DO i = 1, row_length
               rwrk(i,2*(k-1)+1) = uv_comp(i,in_rows,k) *                &
                                      cos_longitude(i,in_rows)
               rwrk(i,2*k)       = uv_comp(i,in_rows,k) *                &
                                      sin_longitude(i,in_rows)
            END DO
         END DO

         CALL global_2d_sums(rwrk, row_length, 1, 0, 0, 2*model_levels, &
                             wrk, proc_row_group)

         DO k = 1, model_levels

            a_np = 2. * wrk(2*(k-1)+1) / global_row_length
            b_np = 2. * wrk(2*k) / global_row_length

            mag_vector_np(k) = SQRT (a_np * a_np + b_np * b_np)

            IF (a_np  ==  0. .and. b_np  ==  0.) THEN
               dir_vector_np(k) = 0.
            ELSE
               dir_vector_np(k) = atan2 (b_np, a_np)
            END IF

         END DO !  k = 1, model_levels

      END IF  !  at_extremity(PNorth)

      IF (at_extremity(PSouth)) THEN

         DO k = 1, model_levels
            DO i = 1, row_length
               rwrk(i,2*(k-1)+1) = uv_comp(i,1,k) * cos_longitude(i,1)
               rwrk(i,2*k)       = uv_comp(i,1,k) * sin_longitude(i,1)
            END DO
         END DO

         CALL global_2d_sums(rwrk, row_length, 1, 0, 0, 2*model_levels, &
                             wrk, proc_row_group)

         DO k = 1, model_levels

            a_sp = 2. * wrk(2*(k-1)+1) / global_row_length
            b_sp = 2. * wrk(2*k) / global_row_length

            mag_vector_sp(k) = SQRT (a_sp * a_sp + b_sp * b_sp)
            IF (a_sp  ==  0. .and. b_sp  ==  0.) THEN
               dir_vector_sp(k) = 0.
            ELSE
               dir_vector_sp(k) = atan2 (b_sp, a_sp)
            END IF

         END DO  !  k = 1, model_levels

      END IF !  at_extremity(PSouth)

! End of routine.
      IF (lhook) CALL dr_hook('POLAR_VECTOR_WIND_N',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Polar_Vector_Wind_n
