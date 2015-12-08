! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine Mpp_tri_solve_setup

      Subroutine Mpp_tri_solve_setup(                                   &
     &                     row_length, rows, model_levels,              &
     &                     j_start, j_end,                              &
     &                     n_proc, n_procx, me,                         &
     &                     proc_row_group,                              &
     &                     a_central_x, a_plus_x, a_minus_x,            &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1,                             &
     &                     bv_soln_1_term2)

! Purpose: Calculates coefficients for MPP periodic tri-diagonal
!          matrix solver.
!
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
!

      Use mpl, Only : &
          mpl_real


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
     &, proc_row_group                                                  &
     &, me                                                              &
     &, j_start, j_end

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
      Real                                                              &
     &  a_central_x(row_length,rows,model_levels)                       &
     &, a_plus_x(row_length+1,rows,model_levels)                        &
     &, a_minus_x(row_length+1,rows,model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  factor_forward(row_length+1,j_start:j_end,model_levels)         &
     &, factor_backward(row_length+1,j_start:j_end,model_levels)        &
     &, bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
     &, bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
     &, recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
     &, bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
     &, bv_factor_forward(2,2*n_procx,rows,model_levels)                &
     &, bv_factor_backward(2*n_procx,rows,model_levels)                 &
     &, recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  bv_soln_n_term(rows, model_levels)                              &
     &, bv_soln_1_term1(rows, model_levels)                             &
     &, bv_soln_1_term2(rows, model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k,d                                                         &
     &, nbv_points                                                      &
     &, ime                                                             &
     &, ibase                                                           &
     &, info                                                            &
     &, dj                                                              &
     &, len

      Real                                                              &
     &  divisor_1

! Local arrays.
      Real                                                              &
     &  bv_a_matrix(-2:2,2*n_procx,rows,model_levels)

! communication buffers, one per x processor
      Real                                                              &
     &  rbuf(6,rows,model_levels,0:n_procx-1)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External routines.

!-----------------------------------------------------------------------
!     Section 1. Forward Gaussian elimination sweep
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('MPP_TRI_SOLVE_SETUP',zhook_in,zhook_handle)
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 2, row_length
! calculate weight to eliminate a_minus_x on this row of the matrix
            factor_forward(i,j,k) = - a_minus_x(i,j,k) /                &
     &                                a_central_x(i-1,j,k)
! update central matrix entry for this row
            a_central_x(i,j,k) = a_central_x(i,j,k) +                   &
     &                           factor_forward(i,j,k) *                &
     &                           a_plus_x(i-1,j,k)
! plus matrix entry is unchanged
! storage space for minus matrix entry is now used to hold value for
! zeroth point
            a_minus_x(i,j,k) = factor_forward(i,j,k) *                  &
     &                         a_minus_x(i-1,j,k)
          End Do
        End Do
      End Do

! set-up reciprocal of a_central_x at this point as a_central_x is
! not changed by any subsequent code. Reciprocal removes divisions
! from code.
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            recip_a_central_x(i,j,k) = 1. / a_central_x(i,j,k)
          End Do
        End Do
      End Do
!-----------------------------------------------------------------------
!     Section 2. Backward Gaussian elimination sweep
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = row_length - 1, 1, -1
! calculate weight to eliminate a_plus_x on this row of the matrix
            factor_backward(i,j,k) = - a_plus_x(i,j,k) *                &
     &                                 recip_a_central_x(i+1,j,k)
! central matrix entry for this row is unchanged as minus entry on
! row below has been removed by forward sweep.
! plus matrix entry is now used to hold value for row_length+1 th
! point.
            a_plus_x(i,j,k) = factor_backward(i,j,k) *                  &
     &                         a_plus_x(i+1,j,k)
! a minus matrix entry is now used to hold value for
! zeroth point and needs updating
            a_minus_x(i,j,k) = a_minus_x(i,j,k) +                       &
     &                         factor_backward(i,j,k) *                 &
     &                         a_minus_x(i+1,j,k)
          End Do
        End Do
      End Do

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
        Do j = 1, rows
          rbuf(1,j,k,ime) = a_minus_x(1,j,k)
          rbuf(2,j,k,ime) = a_central_x(1,j,k)
          rbuf(3,j,k,ime) = a_plus_x(1,j,k)
          rbuf(4,j,k,ime) = a_minus_x(row_length,j,k)
          rbuf(5,j,k,ime) = a_central_x(row_length,j,k)
          rbuf(6,j,k,ime) = a_plus_x(row_length,j,k)
        End Do
      End Do

! broadcast buffer data to all processors on this row.
      len = model_levels * rows * 6


! As GCOM groups are MPI communicators when we are sending data from
! all processors in a group to all processors in the group, we may
! use an MPI collective rather than looping around GCOM broadcasts
      Call mpl_allgather(rbuf(1,1,1,ime), len, MPL_REAL,              &
                         rbuf(1,1,1,  0), len, MPL_REAL,              &
                         proc_row_group,  info)








! initialize bv_a_matrix to zero
      bv_a_matrix(:,:,:,:)=0.0
      bv_a_matrix_0(:,:,:)=0.0
      bv_a_matrix_np1(:,:,:)=0.0

! Copy all boundary data values into new matrix
! limit plus and minus points to range 1 to nbv_points inclusive

      If (n_procx  ==  1) Then
        dj=1
        Do k = 1, model_levels
          Do j = j_start, j_end
! left hand boundary point on processor
            bv_a_matrix(0,dj,j,k) = rbuf(2,j,k,ime)
            bv_a_matrix_np1(dj,j,k) = rbuf(3,j,k,ime)
            bv_a_matrix_0(dj,j,k) = rbuf(1,j,k,ime)
! right hand boundary point on processor
            bv_a_matrix(0,dj+1,j,k) = rbuf(5,j,k,ime)
            bv_a_matrix_np1(dj+1,j,k) = rbuf(6,j,k,ime)
            bv_a_matrix_0(dj+1,j,k) = rbuf(4,j,k,ime)
          End Do
        End Do

      Else
! data from row processor zero goes at start of bv matrix
        dj=1
          Do k = 1, model_levels
            Do j = j_start, j_end
! left hand boundary point on processor
              bv_a_matrix(0,dj,j,k) = rbuf(2,j,k,0)
              bv_a_matrix(2,dj,j,k) = rbuf(3,j,k,0)
              bv_a_matrix_0(dj,j,k) = rbuf(1,j,k,0)
! right hand boundary point on processor
              bv_a_matrix(0,dj+1,j,k) = rbuf(5,j,k,0)
              bv_a_matrix(1,dj+1,j,k) = rbuf(6,j,k,0)
              bv_a_matrix_0(dj+1,j,k) = rbuf(4,j,k,0)
            End Do
          End Do
        dj = (n_procx-1) * 2 + 1
          Do k = 1, model_levels
            Do j = j_start, j_end
! left hand boundary point on processor
              bv_a_matrix(0,dj,j,k) = rbuf(2,j,k,n_procx-1)
              bv_a_matrix_np1(dj,j,k) = rbuf(3,j,k,n_procx-1)
              bv_a_matrix(-1,dj,j,k) = rbuf(1,j,k,n_procx-1)
! right hand boundary point on processor
              bv_a_matrix(0,dj+1,j,k) = rbuf(5,j,k,n_procx-1)
              bv_a_matrix_np1(dj+1,j,k) = rbuf(6,j,k,n_procx-1)
              bv_a_matrix(-2,dj+1,j,k) = rbuf(4,j,k,n_procx-1)
            End Do
          End Do
! data from mid row processors go in the middle
        Do d = 1, n_procx-2
          dj = d * 2 + 1
          Do k = 1, model_levels
            Do j = j_start, j_end
! left hand boundary point on processor
              bv_a_matrix(0,dj,j,k) = rbuf(2,j,k,d)
              bv_a_matrix(2,dj,j,k) = rbuf(3,j,k,d)
              bv_a_matrix(-1,dj,j,k) = rbuf(1,j,k,d)
! right hand boundary point on processor
              bv_a_matrix(0,dj+1,j,k) = rbuf(5,j,k,d)
              bv_a_matrix(1,dj+1,j,k) = rbuf(6,j,k,d)
              bv_a_matrix(-2,dj+1,j,k) = rbuf(4,j,k,d)
            End Do
          End Do
        End Do
      End If ! on n_procx  >   1

      If (n_procx  >   1) Then
! for one processor solution is for a 2x2 matrix and is done entirely
! inside exec part of scheme
!-----------------------------------------------------------------------
!     Section 3.1. Forward Gaussian elimination sweep on bv matrix
!-----------------------------------------------------------------------

        Do k = 1, model_levels
          Do i = j_start, j_end
          Do j = 2, nbv_points - 2
            If (mod(j,2)  ==  0) Then
! elimination from data on even rows
! eliminate sub diagonal entry on row below
              bv_factor_forward(1,j,i,k) = - bv_a_matrix(-1,j+1,i,k) /  &
     &                                     bv_a_matrix(0,j,i,k)
! set entry to zero
              bv_a_matrix(-1,j+1,i,k) = 0.0
! update other entries
! 1. diagonal
              bv_a_matrix(0,j+1,i,k) = bv_a_matrix(0,j+1,i,k)           &
     &                               + bv_factor_forward(1,j,i,k)       &
     &                               * bv_a_matrix(1,j,i,k)
! 2. create 0 point entry
              bv_a_matrix_0(j+1,i,k) = bv_factor_forward(1,j,i,k)       &
     &                               * bv_a_matrix_0(j,i,k)
! eliminate sub-sub-diagonal entry on row 2 below
              bv_factor_forward(2,j,i,k) = - bv_a_matrix(-2,j+2,i,k) /  &
     &                                     bv_a_matrix(0,j,i,k)
! set entry to zero
!              bv_a_matrix(-2,j+2,i,k) = 0.0
! update other entries
! 1. sub-diagonal (is zero to begin with)
              bv_a_matrix(-1,j+2,i,k) = bv_factor_forward(2,j,i,k)      &
     &                                * bv_a_matrix(1,j,i,k)
! 2. create 0 point entry
              bv_a_matrix_0(j+2,i,k) = bv_factor_forward(2,j,i,k)       &
     &                               * bv_a_matrix_0(j,i,k)

            Else
! elimination from data on odd rows
! eliminate sub diagonal entry on row below
              bv_factor_forward(1,j,i,k) = - bv_a_matrix(-1,j+1,i,k) /  &
     &                                     bv_a_matrix(0,j,i,k)
! set entry to zero
!              bv_a_matrix(-1,j+1,i,k) = 0.0
! update other entries
! no sup diagonal on row j so no need to update diagonal entry on j+1
! 1. sup-diagonal entry
              bv_a_matrix(1,j+1,i,k) = bv_a_matrix(1,j+1,i,k)           &
     &                               + bv_factor_forward(1,j,i,k)       &
     &                               * bv_a_matrix(2,j,i,k)
! 2. 0 point entry updating
              bv_a_matrix_0(j+1,i,k) = bv_a_matrix_0(j+1,i,k)           &
     &                               + bv_factor_forward(1,j,i,k)       &
     &                               * bv_a_matrix_0(j,i,k)

            End If
          End Do
          End Do
        End Do
! special case for point nbv_points-1
        Do k = 1, model_levels
          Do i = j_start, j_end
          j = nbv_points - 1
! elimination from data on odd rows
! eliminate sub diagonal entry on row below
            bv_factor_forward(1,j,i,k) = - bv_a_matrix(-1,j+1,i,k) /    &
     &                                     bv_a_matrix(0,j,i,k)
! set entry to zero
!            bv_a_matrix(-1,j+1,i,k) = 0.0
! update other entries
! no sup diagonal on row j so no need to update diagonal entry on j+1
! 1. sup-diagonal entry
            bv_a_matrix_np1(j+1,i,k) = bv_a_matrix_np1(j+1,i,k)         &
     &                               + bv_factor_forward(1,j,i,k)       &
     &                               * bv_a_matrix_np1(j,i,k)
! 2. 0 point entry updating
            bv_a_matrix_0(j+1,i,k) = bv_a_matrix_0(j+1,i,k)             &
     &                               + bv_factor_forward(1,j,i,k)       &
     &                               * bv_a_matrix_0(j,i,k)

          End Do
        End Do

! store recip_bv_a_matrix
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, nbv_points
            recip_bv_a_matrix_diag(i,j,k) = 1./ bv_a_matrix(0,i,j,k)
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.2. Backward Gaussian elimination sweep on bv matrix
!-----------------------------------------------------------------------

        Do k = 1, model_levels
          Do i = j_start, j_end
          Do j = nbv_points - 1, 3, -2
! eliminate sup-sup-diagonal entries from odd rows of matrix
            bv_factor_backward(j,i,k) = - bv_a_matrix(2,j-2,i,k) *      &
     &                                    recip_bv_a_matrix_diag(j,i,k)
! set entry to zero
!            bv_a_matrix(2,j-2,i,k) = 0.0
! update other entries
! 1. np1 point entry creation
            bv_a_matrix_np1(j-2,i,k) = bv_factor_backward(j,i,k)        &
     &                               * bv_a_matrix_np1(j,i,k)
! 2. 0 point entry updating
            bv_a_matrix_0(j-2,i,k) = bv_a_matrix_0(j-2,i,k)             &
     &                               + bv_factor_backward(j,i,k)        &
     &                               * bv_a_matrix_0(j,i,k)
          End Do
          End Do
        End Do

      End If ! on n_procx  >   1

!-----------------------------------------------------------------------
!     Section 3.3. Calculate coefficients for end point solutions.
!-----------------------------------------------------------------------

      Do k = 1, model_levels
      Do j = j_start, j_end
        divisor_1 = 1./ (bv_a_matrix(0,1,j,k) +                         &
     &                     bv_a_matrix_np1(1,j,k))
        bv_soln_n_term(j,k) = 1./ (bv_a_matrix(0,nbv_points,j,k) +      &
     &                          bv_a_matrix_0(nbv_points,j,k))
        bv_soln_1_term1(j,k) = bv_a_matrix_0(1,j,k)                     &
     &                       * bv_soln_n_term(j,k)
        bv_soln_1_term2(j,k) = divisor_1 /                              &
     &                     ( 1. - (bv_a_matrix_0(1,j,k) *               &
     &                             bv_a_matrix_np1(nbv_points,j,k))     &
     &                             * bv_soln_n_term(j,k) * divisor_1)
      End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.4. Store sup-diagonal for
!                  passing out of routine.
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, nbv_points
            bv_a_matrix_sup(i,j,k) = bv_a_matrix(1,i,j,k)
          End Do
        End Do
      End Do

!     end of routine mpp_tri_solve_setup

      IF (lhook) CALL dr_hook('MPP_TRI_SOLVE_SETUP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Mpp_tri_solve_setup
