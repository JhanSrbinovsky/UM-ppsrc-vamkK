! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Calculates pre-conditioning operator applied to field.
!          Set-up of matrix pre-computed.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver

      SUBROUTINE GCR_precon_1_exec(                                     &
     &                     RHS,model_domain,                            &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     a0, a1, factor, Soln,                        &
     &                     offx, offy, at_extremity                     &
     &                     )



      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain

!  Parallel variables
      LOGICAL                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      INTEGER                                                           &
     & offx,offy

      REAL                                                              &
     &  RHS(row_length,rows,model_levels)                               &
                                          ! RIGHT-HAND-SIDE OF EQUATION.
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      REAL                                                              &
     &  a0(row_length,rows,model_levels)                                &
     &, a1(row_length,rows,model_levels)                                &
     &, factor(row_length,rows,model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      REAL                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                ! SOLUTION.

! Local Variables.

      INTEGER                                                           &
     &  i,j,k                                                           &
! cjj bdy.  Additions.
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop

      INTEGER :: jj
      INTEGER :: j1
      INTEGER :: j2
      INTEGER :: len1
      INTEGER :: len2

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Local arrays.

!    No External routines.

!-----------------------------------------------------------------------
!     Section 1.
!-----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('GCR_PRECON_1_EXEC',zhook_in,zhook_handle)
      IF (model_domain  ==  mt_global .OR.                              &
     &    model_domain  ==  mt_bi_cyclic_LAM) THEN
! Solve over full domain
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
      ELSE IF(model_domain  ==  mt_lam)THEN
! Solve over interior points only.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_stop = rows-1
        IF(at_extremity(PEast)) i_stop = row_length-1
        IF(at_extremity(PWest)) i_start = 2
      ELSEIF (model_domain  ==  mt_cyclic_LAM) THEN
! Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_stop = rows-1
      END IF

! Cache-blocked by 4 with some OpenMP added.
      k=1
      DO j = 1, rows
        DO i = 1, row_length
          Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
        END DO
      END DO
      DO k = 2, model_levels
        DO j = 1, rows
          DO i = 1, i_start-1
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          END DO
          DO i = i_stop+1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          END DO
        END DO
        DO j = 1, j_start-1
          DO i = 1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          END DO
        END DO
        DO j = j_stop+1,rows
          DO i = 1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          END DO
        END DO
      END DO

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

      len1=4
      len2=(j_stop-j_start+1)/len1+1
!
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(j,k,i,jj,j1,j2)&
!$OMP& SHARED(RHS,FV_sec_theta_latitude,j_start,j_stop,i_start,i_stop,  &
!$OMP&        a0,Soln,factor,a1,model_levels,len1,len2)
      DO jj=1, len1
        j1=j_start+(jj-1)*len2
        j2=MIN(j_start+jj*len2-1, j_stop)
        DO k= 2, model_levels
          DO j = j1, j2
            DO i = i_start, i_stop
!  the change below causes the results to change
              Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j) -   &
!             Soln(i,j,k) =  Soln(i,j,k) -                              &
     &                          factor(i,j,k)*Soln(i,j,k-1)
            END DO
          END DO
        END DO
        DO j = j1, j2
          DO i = i_start, i_stop
            Soln(i,j,model_levels)  = a0(i,j,model_levels) *            &
     &                                Soln(i,j,model_levels)
          END DO
        END DO
        DO k= model_levels-1, 1, -1
          DO j = j1, j2
            DO i = i_start, i_stop
              Soln(i,j,k)  = a0(i,j,k) *                                &
     &      ( Soln(i,j,k) - a1(i,j,k) * Soln(i,j,k+1) )
            END DO
          END DO
        END DO
      END DO  ! jj
!$OMP END PARALLEL DO

      IF (lhook) CALL dr_hook('GCR_PRECON_1_EXEC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_1_exec
