! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!! Subroutine h_Cubic_v_quintic
      Subroutine h_Cubic_v_quintic                                      &
     &                         (Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j, number_of_inputs,       &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, k_out,                    &
     &                          coeff_z,                                &
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
     &, number_of_inputs ! number of fields to interpolate

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,                             &
                                                   !data to be
     &           -1:dim_k_in+2, number_of_inputs)                       &
                                                   ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number
                                                      ! between 0 & 1

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
     &  i, j, k, index, n ! Loop indices

      Real                                                              &
     &  val_i_minus                                                     &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
     &, coeff_i_minus(dim_i_out, dim_j_out)                             &
     &, coeff_i_zero(dim_i_out, dim_j_out)                              &
     &, coeff_i_plus(dim_i_out, dim_j_out)                              &
     &, coeff_i_plus2(dim_i_out, dim_j_out)                             &
     &, coeff_j_minus(dim_i_out, dim_j_out)                             &
     &, coeff_j_zero(dim_i_out, dim_j_out)                              &
     &, coeff_j_plus(dim_i_out, dim_j_out)                              & 
     &, coeff_j_plus2(dim_i_out, dim_j_out)


      Real :: col_data_i_j_n (-2:3)




                                            ! Horizontally interpolated
                                            ! data ready for vertical
                                            ! interpolation.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None.

      IF (lhook) CALL dr_hook('H_CUBIC_V_QUINTIC',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 0.   Calculate all interpolation weights
! ----------------------------------------------------------------------


!$OMP PARALLEL DO SCHEDULE(STATIC,1) DEFAULT(NONE)                      &
!$OMP& PRIVATE(i,j,k,index,n,val_i_minus,val_i,val_i_plus,val_i_plus2,  &
!$OMP& phi_cubed,phi_sq,lambda_cubed,lambda_sq,col_data_i_j_n,          &
!$OMP& coeff_j_minus, coeff_j_plus2,coeff_j_plus,coeff_j_zero,          &
!$OMP& coeff_i_plus2, coeff_i_plus,coeff_i_minus,coeff_i_zero)          &
!$OMP& SHARED(weight_phi,weight_lambda, dim_k_out,dim_j_out,dim_i_out,  &
!$OMP& number_of_inputs,Ext_Data,i_out,j_out,k_out,coeff_z,Data_out)

      Do k = 1, dim_k_out
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
!             Calculate interpolation weights in the j direction
            phi_sq    = weight_phi (i,j,k) * weight_phi (i,j,k)
            phi_cubed = weight_phi (i,j,k) * phi_sq

            coeff_j_plus2(i,j) = phi_cubed - weight_phi(i,j,k)
            coeff_j_plus(i,j)  = 0.5 *                                  &
     &                      (coeff_j_plus2(i,j) - phi_sq -              &
     &                       weight_phi(i,j,k))
            coeff_j_zero(i,j)  = 0.5 *                                  &
     &           (coeff_j_plus2(i,j) - 2. * phi_sq  + 2.)
            coeff_j_plus2(i,j) = one_sixth * coeff_j_plus2(i,j)
            coeff_j_minus(i,j) = one_sixth *                            &
     &                (phi_cubed - 3.*phi_sq + 2. * weight_phi(i,j,k))
            lambda_sq = weight_lambda(i,j,k) * weight_lambda(i,j,k)
            lambda_cubed = weight_lambda (i,j,k) * lambda_sq

            coeff_i_plus2(i,j)  = lambda_cubed -                        &
     &               weight_lambda(i,j,k)
            coeff_i_plus(i,j)   = 0.5 *                                 &
     &             (coeff_i_plus2(i,j) - lambda_sq -                    &
     &                   weight_lambda(i,j,k))
            coeff_i_zero(i,j)   = 0.5 *                                 &
     &             (coeff_i_plus2(i,j) - 2.*lambda_sq + 2.)
            coeff_i_plus2(i,j)  = one_sixth * coeff_i_plus2(i,j)
            coeff_i_minus(i,j)  = one_sixth *                           &
     &       (lambda_cubed - 3.*lambda_sq + 2. * weight_lambda(i,j,k))

          End Do ! i 
        End Do ! j
! ----------------------------------------------------------------------
! Section 1.   Perform cubic Lagrange Interpolation in j direction.
! ----------------------------------------------------------------------

        Do n = 1, number_of_inputs
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out  

!         loop over 6 levels required for quintic vertical interpolation
              Do index = -2, 3

!             Interpolate data to the 4 i points needed to do the
!             i direction interpolation
                val_i_minus = coeff_j_plus2(i,j) *                      &
     &    Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,k_out(i,j,k)+index,n) &
     &                       - coeff_j_plus(i,j) *                      &
     &    Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,k_out(i,j,k)+index,n) &
     &                       + coeff_j_zero(i,j) *                      &
     &    Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),  k_out(i,j,k)+index,n) &
     &                         - coeff_j_minus(i,j) *                   &
     &    Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,k_out(i,j,k)+index,n)

                val_i = coeff_j_plus2(i,j) *                            &
     &    Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,k_out(i,j,k)+index,n)   &
     &                - coeff_j_plus(i,j) *                             &
     &    Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+index,n)   &
     &                + coeff_j_zero(i,j) *                             &
     &    Ext_Data (i_out(i,j,k),j_out(i,j,k)  ,k_out(i,j,k)+index,n)   &
     &                - coeff_j_minus(i,j) *                            &
     &    Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,k_out(i,j,k)+index,n)

                val_i_plus = coeff_j_plus2(i,j) *                       &
     &    Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,k_out(i,j,k)+index,n) &
     &                     - coeff_j_plus(i,j) *                        &
     &    Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+index,n) &
     &                     + coeff_j_zero(i,j) *                        &
     &    Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)  ,k_out(i,j,k)+index,n) &
     &                     - coeff_j_minus(i,j) *                       &
     &    Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,k_out(i,j,k)+index,n)

                val_i_plus2 = coeff_j_plus2(i,j) *                      &
     &    Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,k_out(i,j,k)+index,n) &
     &                      - coeff_j_plus(i,j) *                       &
     &    Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,k_out(i,j,k)+index,n) &
     &                      + coeff_j_zero(i,j) *                       &
     &    Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)  ,k_out(i,j,k)+index,n) &
     &                      - coeff_j_minus(i,j) *                      &
     &    Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,k_out(i,j,k)+index,n)

! ----------------------------------------------------------------------
! Section 2.   Perform cubic Lagrange Interpolation in i direction.
! ----------------------------------------------------------------------

! Interpolate data to get the data at the point for the k interpolation

                col_data_i_j_n (index) =                                &
     &                             coeff_i_plus2(i,j) * val_i_plus2     &
     &                           - coeff_i_plus(i,j)  * val_i_plus      &
     &                           + coeff_i_zero(i,j)  * val_i           &
     &                           - coeff_i_minus(i,j) * val_i_minus
              End Do ! end loop over index

! ----------------------------------------------------------------------
! Section 3.   Perform quintic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

              Data_out(i,j,k,n) = coeff_z(i,j,k,-2) *col_data_i_j_n (-2)&
     &                          + coeff_z(i,j,k,-1) *col_data_i_j_n (-1)&
     &                          + coeff_z(i,j,k,0)  *col_data_i_j_n (0) &
     &                          + coeff_z(i,j,k,1)  *col_data_i_j_n (1) &
     &                          + coeff_z(i,j,k,2)  *col_data_i_j_n (2) &
     &                          + coeff_z(i,j,k,3)  *col_data_i_j_n (3)
            End Do  ! i 
          End Do  ! j
        End Do  ! n = 1, number_of_inputs
      End Do  ! k
!$OMP END PARALLEL DO

!     End of routine.
      IF (lhook) CALL dr_hook('H_CUBIC_V_QUINTIC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE h_Cubic_v_quintic
