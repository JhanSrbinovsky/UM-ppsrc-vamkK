! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine Mpp_tri_solve_exec

      Subroutine Mpp_tri_solve_exec(                                    &
     &                     row_length, rows, model_levels,              &
     &                     off_x, off_y, j_start, j_end,                &
     &                     n_proc, n_procx, me,                         &
     &                     proc_row_group,                              &
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1,                             &
     &                     bv_soln_1_term2,                             &
     &                     rhs_xig, sss, soln)

! Purpose:
!          Calculates ADI pre-conditioning matrix coefficients.
!          This code only for global model at the moment.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      Use mpl, Only : &
          MPL_REAL


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit none

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, me                                                              &
     &, proc_row_group                                                  &
     &, off_x                                                           &
     &, off_y                                                           &
     &, j_start, j_end

      real con1,con2
      Real                                                              &
     &  a_plus_x(row_length+1,rows,model_levels)                        &
     &, a_minus_x(row_length+1,rows,model_levels)                       &
     &, rhs_xig(row_length,rows,model_levels)

      real,intent(in),                                                  &
     &     dimension(row_length+1,j_start:j_end,model_levels) :: sss

      Real                                                              &
     &  factor_forward(row_length+1,j_start:j_end,model_levels)         &
     &, factor_backward(row_length+1,j_start:j_end,model_levels)        &
     &, bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
     &, bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
     &, bv_factor_forward(2,2*n_procx,rows,model_levels)                &
     &, bv_factor_backward(2*n_procx,rows,model_levels)                 &
     &, recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
     &, bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
     &, recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  bv_soln_n_term(rows, model_levels)                              &
     &, bv_soln_1_term1(rows, model_levels)                             &
     &, bv_soln_1_term2(rows, model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  soln(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k,d,jj                                                      &
     &, nbv_points                                                      &
     &, ime                                                             &
     &, ibase                                                           &
     &, info                                                            &
     &, len                                                             &
     &, start_point

      Real :: rhs_x_1    ! unrolling temporary
      Real :: rhs_x_2    ! unrolling temporary
      Real :: rhs_x_3    ! unrolling temporary
      Real :: rhs_x_4    ! unrolling temporary

! Local arrays.
      Real                                                              &
     &  bv_rhs_x(2*n_procx,rows,model_levels)                           &
     &, bv_soln(0:2*n_procx+1,rows,model_levels)                        &
     &, rbuf(2,rows,model_levels,0:n_procx-1)

      Real rhs_xx(row_length+1,j_start:j_end,model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External routines.
!-----------------------------------------------------------------------
!     Section 1. Forward Gaussian elimination sweep
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('MPP_TRI_SOLVE_EXEC',zhook_in,zhook_handle)

! The following code snippet is the basic algorithm that has been optimised
! for vector and scalar platforms below and is given for reference.
!      Do k = 1, model_levels    
!        Do j = j_start, j_end    
!          rhs_xx(1,j,k)=sss(1,j,k)    
!        End Do    
!     
!        Do j = j_start, j_end    
!          Do i = 2, row_length    
!! update right-hand-side term    
!            rhs_xx(i,j,k) = sss(i,j,k) + factor_forward(i,j,k) *        &    
!                                          rhs_xx(i-1,j,k)    
!          End Do    
!    
!          Do i = row_length - 1, 1, -1    
!! update right-hand-side term    
!            rhs_xx(i,j,k) = rhs_xx(i,j,k) + factor_backward(i,j,k) *    &    
!                                          rhs_xx(i+1,j,k)    
!          End Do    
!        End Do    
!      End Do 

! Hand unrolling by 4 with some
! OpenMP added. Using the UNROLL directive doesn't seem to help unfortunately
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(i, j, k, jj, rhs_x_1, rhs_x_2, rhs_x_3, rhs_x_4)         &
!$OMP& SHARED(model_levels, j_start, j_end, row_length, rhs_xx, sss,    &
!$OMP&        factor_forward, factor_backward)
      Do k = 1, model_levels
        Do j = j_start, j_end
          rhs_xx(1,j,k)=sss(1,j,k)
        End Do
 
        Do j = j_start, j_end-3, 4
! The rhs_x_(1-4) variables are the 4-strided rows of the first column
! of rhs_xx at this point
          rhs_x_1 = rhs_xx(1,j  ,k)
          rhs_x_2 = rhs_xx(1,j+1,k)
          rhs_x_3 = rhs_xx(1,j+2,k)
          rhs_x_4 = rhs_xx(1,j+3,k)

          Do i = 2, row_length
! update right-hand-side term
            rhs_x_1 = sss(i,j  ,k) + factor_forward(i,j  ,k) * rhs_x_1
            rhs_x_2 = sss(i,j+1,k) + factor_forward(i,j+1,k) * rhs_x_2
            rhs_x_3 = sss(i,j+2,k) + factor_forward(i,j+2,k) * rhs_x_3
            rhs_x_4 = sss(i,j+3,k) + factor_forward(i,j+3,k) * rhs_x_4
            rhs_xx(i,j  ,k) = rhs_x_1
            rhs_xx(i,j+1,k) = rhs_x_2
            rhs_xx(i,j+2,k) = rhs_x_3
            rhs_xx(i,j+3,k) = rhs_x_4
          End Do
 
! Now start the backward sweep - so start with data at end of rows
! 4-strided again.
          rhs_x_1 = rhs_xx(row_length,j  ,k)
          rhs_x_2 = rhs_xx(row_length,j+1,k)
          rhs_x_3 = rhs_xx(row_length,j+2,k)
          rhs_x_4 = rhs_xx(row_length,j+3,k)

          Do i = row_length - 1, 1, -1
! update right-hand-side term
            rhs_x_1  = rhs_xx(i,j  ,k) + factor_backward(i,j  ,k) * &
     &                                          rhs_x_1
            rhs_x_2 =  rhs_xx(i,j+1,k) + factor_backward(i,j+1,k) * &
     &                                          rhs_x_2
            rhs_x_3  = rhs_xx(i,j+2,k) + factor_backward(i,j+2,k) * &
     &                                          rhs_x_3
            rhs_x_4 =  rhs_xx(i,j+3,k) + factor_backward(i,j+3,k) * &
     &                                          rhs_x_4
            rhs_xx(i,j  ,k) = rhs_x_1
            rhs_xx(i,j+1,k) = rhs_x_2
            rhs_xx(i,j+2,k) = rhs_x_3
            rhs_xx(i,j+3,k) = rhs_x_4
          End Do

        jj = j
        End Do

! We've done the main calculation of rhs_xx, but as we've strided in
! blocks of 4 there may be a few j elements not yet covered. These
! will be mopped up starting after the finishing value in the previous
! loop (stored in jj). The next value we need to compute is at jj+4
! as the previous values have already been calculated in the strided loop.
        jj = jj + 4
        Do j=jj, j_end
          Do i = 2, row_length
            rhs_xx(i,j  ,k) = sss(i,j  ,k) + factor_forward(i,j  ,k) * &
     &                                       rhs_xx(i-1,j  ,k)
          End Do

          Do i = row_length - 1, 1, -1
            rhs_xx(i  ,j,k) = rhs_xx(i  ,j,k) + factor_backward(i  ,j,k) * &
     &                                          rhs_xx(i+1,j  ,k)
          End Do
        End Do
      End Do
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!     Section 3. Solve for boundary points
!-----------------------------------------------------------------------
! calculate position on row of processors, 0 = start, n_procx-1 = end
      ibase = (me/n_procx) * n_procx
      ime = me - ibase

! number of boundary data points is twice number of processors
      nbv_points = 2 * n_procx

! copy data to be broadcast from this processor into the buffer
      Do k = 1, model_levels
        Do j = j_start, j_end
          rbuf(1,j,k,ime) = rhs_xx(1,j,k)
          rbuf(2,j,k,ime) = rhs_xx(row_length,j,k)
        End Do
      End Do

! broadcast buffer data to all processors on this row.
      len = model_levels * rows * 2

      if(n_procx > 1)  then
! As GCOM groups are MPI communicators when we are sending data from
! all processors in a group to all processors in the group, we may
! use an MPI collective rather than looping around GCOM broadcasts
        Call mpl_allgather(rbuf(1,1,1,ime), len, MPL_REAL,              &
                           rbuf(1,1,1,  0), len, MPL_REAL,              &
                           proc_row_group,  info)
      end if

! Copy all boundary data values into new matrix
      Do d = 0, n_procx-1
        i = 2 * d
        Do k = 1, model_levels
          Do j = j_start, j_end
            bv_rhs_x(i+1,j,k) = rbuf(1,j,k,d)
            bv_rhs_x(i+2,j,k) = rbuf(2,j,k,d)
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.1. Forward Gaussian elimination sweep on bv rhs
!-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& PRIVATE(i, j, k, start_point, con1, con2)                        &
!$OMP& SHARED(j_start, j_end, nbv_points, model_levels, bv_rhs_x,       &
!$OMP&        bv_factor_forward, bv_factor_backward, bv_soln,           &
!$OMP&        bv_soln_1_term1, bv_soln_1_term2, bv_a_matrix_np1,        &
!$OMP&        bv_soln_n_term, bv_a_matrix_0, recip_bv_a_matrix_diag,    &
!$OMP&        soln, rhs_xx, a_minus_x, a_plus_x, recip_a_central_x,     &
!$OMP&        row_length, ime, bv_a_matrix_sup)

! Note that in this parallel section we have many OMP Do loops with 
! NOWAIT clauses. The NOWAITs are only valid as we have a STATIC
! SCHEDULE and the parallelised loops are over the same ranges in
! all cases. This ensures a thread will get the same ranges of levels
! throughout and dependencies between loops are thus satisfied.
!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
        Do j = j_start, j_end
        Do i = 2, nbv_points - 2, 2
! elimination from data on even rows
          bv_rhs_x(i+1,j,k) = bv_rhs_x(i+1,j,k)                         &
     &                      + bv_factor_forward(1,i,j,k)                &
     &                      * bv_rhs_x(i,j,k)
          bv_rhs_x(i+2,j,k) = bv_rhs_x(i+2,j,k)                         &
     &                      + bv_factor_forward(2,i,j,k)                &
     &                      * bv_rhs_x(i,j,k)

! elimination from data on odd rows
          bv_rhs_x(i+2,j,k) = bv_rhs_x(i+2,j,k)                         &
     &                      + bv_factor_forward(1,i+1,j,k)              &
     &                      * bv_rhs_x(i+1,j,k)
        End Do
        End Do
      End Do
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!     Section 3.2. Backward Gaussian elimination sweep on bv rhs
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
      Do j = j_start, j_end
        Do i = nbv_points - 1, 3, -2
           bv_rhs_x(i-2,j,k) = bv_rhs_x(i-2,j,k)                        &
     &                         + bv_factor_backward(i,j,k)              &
     &                         * bv_rhs_x(i,j,k)
        End Do
      End Do
      End Do
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!     Section 3.3. Solution for end points 1 and nbv_points
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
      Do j = j_start, j_end
        bv_soln(1,j,k) = ( bv_rhs_x(1,j,k)                              &
     &                     - bv_soln_1_term1(j,k) *                     &
     &                       bv_rhs_x(nbv_points,j,k) )                 &
     &                     * bv_soln_1_term2(j,k)
        bv_soln(nbv_points,j,k) = ( bv_rhs_x(nbv_points,j,k) -          &
     &                                bv_a_matrix_np1(nbv_points,j,k) * &
     &                                bv_soln(1,j,k) ) *                &
     &                              bv_soln_n_term(j,k)
      End Do
      End Do
!$OMP END DO NOWAIT

! set 0, and nbv_points+1 solutions, ie make solution periodic
!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
        Do j = j_start, j_end
          bv_soln(0,j,k) = bv_soln(nbv_points,j,k)
          bv_soln(nbv_points+1,j,k) = bv_soln(1,j,k)
        End Do
      End Do
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!     Section 3.4. Find solution in inner of matrix
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
        Do j = j_start, j_end
! first all odd rows
        Do i = 3, nbv_points - 1, 2
          bv_soln(i,j,k) = ( bv_rhs_x(i,j,k)                            &
     &                     -   bv_a_matrix_0(i,j,k) * bv_soln(0,j,k)    &
     &                     -   bv_a_matrix_np1(i,j,k) *                 &
     &                         bv_soln(nbv_points+1,j,k) )              &
     &                     * recip_bv_a_matrix_diag(i,j,k)
        End Do
! now all even rows
        Do i = 2, nbv_points - 2, 2
          bv_soln(i,j,k) = ( bv_rhs_x(i,j,k)                            &
     &                     -   bv_a_matrix_sup(i,j,k) * bv_soln(i+1,j,k)&
     &                     -   bv_a_matrix_0(i,j,k) * bv_soln(0,j,k)    &
     &                     -   bv_a_matrix_np1(i,j,k) *                 &
     &                         bv_soln(nbv_points+1,j,k) )              &
     &                     * recip_bv_a_matrix_diag(i,j,k)
        End Do
        End Do
      End Do
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!     Section 3.5.  Copy bv_soln into main solution
!-----------------------------------------------------------------------

      start_point = 2 * ime
!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
        Do j = j_start, j_end
          soln(0,j,k) = bv_soln(start_point,j,k)
          soln(1,j,k) = bv_soln(start_point+1,j,k)
          soln(row_length,j,k) = bv_soln(start_point+2,j,k)
          soln(row_length+1,j,k) = bv_soln(start_point+3,j,k)
        End Do
      End Do
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!     Section 4. Solve for main solution not at boundaries
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels
        Do j = j_start, j_end
        con1=soln(0,j,k)
        con2=soln(row_length+1,j,k)
!dir$ unroll
          Do i = 2, row_length-1
            soln(i,j,k) = ( rhs_xx(i,j,k)                               &
     &                         - a_minus_x(i,j,k) * con1                &
     &                         - a_plus_x(i,j,k) * con2 )               &
     &                     * recip_a_central_x(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO
!$OMP END PARALLEL

!     end of routine mpp_tri_solve_exec
      IF (lhook) CALL dr_hook('MPP_TRI_SOLVE_EXEC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Mpp_tri_solve_exec
