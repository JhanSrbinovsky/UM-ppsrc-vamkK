! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Cubic_Lagrange
!

      Subroutine Cubic_Lagrange(                                        &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j, number_of_inputs,       &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, k_out,                    &
     &                          row_length_in, rows_in,                 &
     &                          lambda_rm, lambda_rp, phi_rm, phi_rp,   &
     &                          recip_lambda_m, recip_lambda_0,         &
     &                          recip_lambda_p, recip_lambda_p2,        &
     &                          recip_phi_m, recip_phi_0,               &
     &                          recip_phi_p, recip_phi_p2,              &
     &                          coeff_z, L_regular,                     &
     &                          model_domain,                           &
     &                          at_extremity, n_procx, n_procy,         &
     &                          global_row_length, global_rows,         &
     &                          proc_col_group, proc_row_group,         &
     &                          datastart,                              &
     &                          Data_out)

! Purpose:
!          Performs cubic Lagrange interpolation of the input field to a
!          set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out.
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


      USE um_types, ONLY: integer32
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_regular                                                       &
                         ! false if variable resolution
     &, at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, number_of_inputs                                                &
                            ! number of fields to interpolate
     &, row_length_in                                                   &
                            ! for lambda dynamic arrays
     &, rows_in                                                         &
                            ! for phi dynamic arrays
     &, model_domain                                                    &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, proc_col_group                                                  &
     &, global_row_length                                               &
                            ! number of points in domain row
     &, global_rows                                                     &
                            ! number of rows in domain
     &, datastart(3)

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,                             &
                                                   !data to be
     &           -1:dim_k_in+2, number_of_inputs)                       &
                                                   ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

        Real                                                            &
     &  lambda_rm (1-halo_i : row_length_in+halo_i)                     &
     &, lambda_rp (1-halo_i : row_length_in+halo_i)                     &
     &, phi_rm    ( 1-halo_i : row_length_in + halo_i                   &
     &,             1-halo_j : rows_in + halo_j )                       &
     &, phi_rp    ( 1-halo_i : row_length_in + halo_i                   &
     &,             1-halo_j : rows_in + halo_j )                       &
     &, recip_lambda_m (1-halo_i : row_length_in+halo_i)                &
     &, recip_lambda_0 (1-halo_i : row_length_in+halo_i)                &
     &, recip_lambda_p (1-halo_i : row_length_in+halo_i)                &
     &, recip_lambda_p2(1-halo_i : row_length_in+halo_i)                &
     &, recip_phi_m ( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in + halo_j )                     &
     &, recip_phi_0 ( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in + halo_j )                     &
     &, recip_phi_p ( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in + halo_j )                     &
     &, recip_phi_p2( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in+ halo_j )

      INTEGER (KIND=integer32) ::                                       &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
     &, k_out (dim_i_out, dim_j_out, dim_k_out)       ! lies between it
                                                      ! and it+1

! Arguments with Intent IN/OUT. ie: Output variables.
      Real                                                              &
     &  coeff_z (dim_i_out, dim_j_out, dim_k_out,-2:3)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (dim_i_out, dim_j_out, dim_k_out, number_of_inputs)
                      ! data interpolated to desired locations.

! Local Parameters

      Real one_sixth
      Parameter ( one_sixth = 1./6. )

! Local Variables.

      Integer                                                           &
     &  i, j, k, index, n                                               &
                          ! Loop indices
     &, i_index, j_index


      Real                                                              &
     &  val_i_minus                                                     &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
                            ! x1-x6 are temporary variables
     &, x1, x3, x4, x5, x6                                              &
     &, coeff_i_minus(dim_i_out, dim_j_out)                             &
     &, coeff_i_zero(dim_i_out, dim_j_out)                              &
     &, coeff_i_plus(dim_i_out, dim_j_out)                              &
     &, coeff_i_plus2(dim_i_out, dim_j_out)                             &
     &, coeff_j_minus(dim_i_out, dim_j_out)                             & 
     &, coeff_j_zero(dim_i_out, dim_j_out)                              &
     &, coeff_j_plus(dim_i_out, dim_j_out)                              &
     &, coeff_j_plus2(dim_i_out, dim_j_out)

      Real                                                              &
     &  col_data (dim_i_out, dim_j_out, -1:2, number_of_inputs)

                                             ! Horizontally interpolated
                                             ! data ready for vertical
                                             ! interpolation.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None.

! ----------------------------------------------------------------------
! Section 0.   Calculate Horizontal interpolation weights
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CUBIC_LAGRANGE',zhook_in,zhook_handle)


! Interpolate data to get the data at the point for the k interpolation
!! rest of benchmark changes change results
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED)                    &
!$OMP& PRIVATE(i,j,k,index,n,val_i_plus2,val_i_plus,val_i,val_i_minus,&
!$OMP& col_data, x1,x3,x4,x5,x6,i_index,j_index,coeff_j_plus2,        &
!$OMP& coeff_j_plus,coeff_j_zero,coeff_j_minus, coeff_i_plus2,        &
!$OMP& coeff_i_plus,coeff_i_zero,coeff_i_minus, lambda_sq, phi_cubed, &
!$OMP& lambda_cubed, phi_sq)         

      Do k = 1, dim_k_out
        j = 1

        If ( L_regular ) Then

          Do j = 1, dim_j_out
          Do i = 1, dim_i_out




            phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
            phi_cubed = weight_phi (i,j,k) * phi_sq

            coeff_j_plus2(i,j) = phi_cubed - weight_phi(i,j,k)
            coeff_j_plus(i,j)  =                                        &
              &  0.5 * (coeff_j_plus2(i,j) - phi_sq                     &
              &                     - weight_phi(i,j,k))
            coeff_j_zero(i,j)  = 0.5 * (coeff_j_plus2(i,j)              &
              &                             - 2. * phi_sq  + 2.)
            coeff_j_plus2(i,j) = one_sixth * coeff_j_plus2(i,j)
            coeff_j_minus(i,j) = one_sixth * (phi_cubed                 &
              & - 3.*phi_sq + 2. * weight_phi(i,j,k))

            lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
            lambda_cubed = weight_lambda (i,j,k) * lambda_sq

            coeff_i_plus2(i,j)  = lambda_cubed - weight_lambda(i,j,k)
            coeff_i_plus(i,j)   =                                       &
              & 0.5 * (coeff_i_plus2(i,j) - lambda_sq                   &
              & - weight_lambda(i,j,k))
            coeff_i_zero(i,j)   = 0.5 * (coeff_i_plus2(i,j)             &
              & - 2.*lambda_sq + 2.)
            coeff_i_plus2(i,j)  = one_sixth * coeff_i_plus2(i,j)
            coeff_i_minus(i,j)  = one_sixth * (lambda_cubed             &
              & - 3.*lambda_sq + 2. * weight_lambda(i,j,k))

         End Do  !  i loop

         END DO  !  j loop


       Else ! variable resolution

         Do j = 1, dim_j_out
         Do i = 1, dim_i_out



           i_index = i_out(i,j,k)
           j_index = j_out(i,j,k)

           x1 = weight_phi(i,j,k) + phi_rm(i_index,j_index)
           x3 = weight_phi(i,j,k) - 1.0
           x4 = x3 - phi_rp(i_index,j_index)          
           x5 = x3 * x4
           x6 = x1 * weight_phi(i,j,k)

           coeff_j_plus2(i,j) = x6 * x3 * recip_phi_p2(i_index,j_index)          
           coeff_j_plus(i,j)  = x6 * x4 * recip_phi_p(i_index,j_index)          
           coeff_j_zero(i,j)  = x1 * x5 * recip_phi_0(i_index,j_index)          
           coeff_j_minus(i,j) = weight_phi(i,j,k) * x5 *                &
             & recip_phi_m(i_index,j_index)          

           x1 = weight_lambda(i,j,k) + lambda_rm(i_index)
           x3 = weight_lambda(i,j,k) - 1.0
           x4 = x3 - lambda_rp(i_index)
           x5 = x3 * x4
           x6 = x1 * weight_lambda(i,j,k)

           coeff_i_plus2(i,j) = x6 * x3 * recip_lambda_p2(i_index)      
           coeff_i_plus(i,j)  = x6 * x4 * recip_lambda_p(i_index)     
           coeff_i_zero(i,j)  = x1 * x5 * recip_lambda_0(i_index)     
           coeff_i_minus(i,j) = weight_lambda(i,j,k) * x5 *             &
             & recip_lambda_m(i_index)     

        End Do  !  i loop

        END DO  !  j loop

        
      End If

        Do n = 1, number_of_inputs

          Do j = 1, dim_j_out
          Do i = 1, dim_i_out




!CDIR UNROLL=4
            Do index = -1, 2
                val_i_minus         = coeff_j_plus2(i,j) *              &
     &        Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,                  &
     &                  k_out(i,j,k)+index,n)                           &
     &                            - coeff_j_plus(i,j) *                 &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j) *                    &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j) *                   &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,                     &
     &               k_out(i,j,k)+index,n)

                val_i              = coeff_j_plus2(i,j) *               &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j) *                    &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j) *                    &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),                         &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j) *                   &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,                       &
     &               k_out(i,j,k)+index,n)

                val_i_plus         = coeff_j_plus2(i,j) *               &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j) *                    &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j) *                    &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j) *                   &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,                     &
     &               k_out(i,j,k)+index,n)

                val_i_plus2        = coeff_j_plus2(i,j) *               &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j) *                    &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j) *                    &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j) *                   &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,                     &
     &               k_out(i,j,k)+index,n)

! Interpolate data to the 4 i points needed to do the i direction
! interpolation
! ----------------------------------------------------------------------
! Section 2.   Perform cubic Lagrange Interpolation in i direction.
! ----------------------------------------------------------------------

                col_data (i,j,index,n) = coeff_i_plus2(i,j)             &
     &                                   * val_i_plus2                  &
     &                                 - coeff_i_plus(i,j)              &
     &                                   * val_i_plus                   &
     &                                 + coeff_i_zero(i,j)              &
     &                                   * val_i                        &
     &                                 - coeff_i_minus(i,j)             &
     &                                   * val_i_minus
            End Do ! index = -1, 2
          End Do  !  i loop

          End Do  !  j loop

        End Do      !n

! ----------------------------------------------------------------------
! Section 3.   Perform cubic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Interpolate data
        Do n = 1, number_of_inputs

          Do j = 1, dim_j_out
          Do i = 1, dim_i_out




            Data_out (i,j,k,n) = coeff_z(i,j,k,-1) *                    &
     &                             col_data (i,j,-1,n) +                &
     &                             coeff_z(i,j,k,0) *                   &
     &                             col_data (i,j,0,n) +                 &
     &                             coeff_z(i,j,k,1) *                   &
     &                             col_data (i,j,1,n) +                 &
     &                             coeff_z(i,j,k,2) *                   &
     &                             col_data (i,j,2,n)

          End Do !  i loop

          End Do !  j loop

        End Do   !  n = number_of_inputs

      End Do   !  k = 1, dim_k_out
!$OMP END PARALLEL DO

      IF (lhook) CALL dr_hook('CUBIC_LAGRANGE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Cubic_Lagrange
