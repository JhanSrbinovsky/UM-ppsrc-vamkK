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

      SUBROUTINE GCR_precon_1_exec_2B(                                  &
     &                     RHS,model_domain,                            &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     a0, a1, factor, Soln,                        &
     &                     offx, offy,                                  &
     &                     i_start, i_stop, j_start, j_stop             &
     &                     )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain

!  Parallel variables

      Integer                                                           &
     &  offx, offy                                                      &
     &, i_start, i_stop                                                 &
                                   ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop            ! loop bounds set in PE_Helmholtz

      Real                                                              &
     &  RHS(row_length,rows,model_levels)                               &
                                          ! RIGHT-HAND-SIDE OF EQUATION.
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
     &  a0(row_length,rows,model_levels)                                &
     &, a1(row_length,rows,model_levels)                                &
     &, factor(row_length,rows,model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  Soln(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
                ! SOLUTION.

! Local Variables.

      Integer                                                           &
     &  i,j,k

      Integer :: jj
      Integer :: j1
      Integer :: j2
      Integer :: len1
      Integer :: len2

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Local arrays.

!    No External routines.

!-----------------------------------------------------------------------
!     Section 1.
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_PRECON_1_EXEC_2B',zhook_in,zhook_handle)

! Cache-blocked by 4 with some OpenMP added.
      k=1
      Do j = 1, rows
        Do i = 1, row_length
          Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
        End Do
      End Do
      Do k = 2, model_levels
        Do j = 1, rows
          Do i = 1, i_start-1
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          End Do
          Do i = i_stop+1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          End Do
        End Do
        Do j = 1, j_start-1
          Do i = 1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          End Do
        End Do
        Do j = j_stop+1,rows
          Do i = 1, row_length
            Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j)
          End Do
        End Do
      End Do

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

      len1=4
      len2=(j_stop-j_start+1)/len1+1
!
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP& PRIVATE(j,k,i,jj,j1,j2) SHARED(j_start,j_stop,i_start,i_stop,a0, &
!$OMP& Soln,RHS,FV_sec_theta_latitude,factor,a1,model_levels,len1,len2)
      Do jj=1, len1
        j1=j_start+(jj-1)*len2
        j2=min(j_start+jj*len2-1, j_stop)
        Do k= 2, model_levels
          Do j = j1, j2
            Do i = i_start, i_stop
! the change in the line below changes the results
              Soln(i,j,k) = RHS(i,j,k) * FV_sec_theta_latitude(i,j) -   &
     &                          factor(i,j,k)*Soln(i,j,k-1)
            End Do
          End Do
        End Do
          Do j = j1, j2
          Do i = i_start, i_stop
            Soln(i,j,model_levels)  = a0(i,j,model_levels) *            &
     &                                Soln(i,j,model_levels)
          End Do
        End Do
        Do k= model_levels-1, 1, -1
          Do j = j1, j2
            Do i = i_start, i_stop
              Soln(i,j,k)  = a0(i,j,k) *                                  &
     &      ( Soln(i,j,k) - a1(i,j,k) * Soln(i,j,k+1) )
            End Do
          End Do
        End Do      ! k =model_levels
      End Do      ! jj loop
!$OMP END PARALLEL DO

      IF (lhook) CALL dr_hook('GCR_PRECON_1_EXEC_2B',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_precon_1_exec_2B
