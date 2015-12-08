! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine ECMWF_mono_quasi_cubic

      Subroutine ECMWF_mono_quasi_cubic(                                &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j, number_of_inputs,       &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, k_out,                    &
     &                          coeff_z,                                &
     &                          Data_out)

! Purpose:
!          Performs monotone quasi-cubic interpolation, as developed
!          by Mariano Hortal at ECMWF, of the input field to a
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
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0               Original code
!  5.3     11/07/01  bug fix             Andy Malcolm
!  6.0  18/08/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
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
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1
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

! Local Variables.

      Integer                                                           &
     &  i, j, k, index, n ! Loop indices

      Real                                                              &
     &  one_sixth                                                       &
     &, val_i_minus                                                     &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
     &, coeff_i_minus  (dim_i_out, dim_j_out)                           &
     &, coeff_i_zero   (dim_i_out, dim_j_out)                           &
     &, coeff_i_plus   (dim_i_out, dim_j_out)                           &
     &, coeff_i_plus2  (dim_i_out, dim_j_out)                           &
     &, coeff_minus(dim_i_out,dim_j_out)                                &
     &, coeff_zero(dim_i_out,dim_j_out)                                 &
     &, coeff_plus(dim_i_out,dim_j_out)                                 &
     &, coeff_plus2(dim_i_out,dim_j_out)                                &
     &, coeff_j_minus  (dim_i_out, dim_j_out)                           &
     &, coeff_j_zero   (dim_i_out, dim_j_out)                           &
     &, coeff_j_plus   (dim_i_out, dim_j_out)                           &
     &, coeff_j_plus2  (dim_i_out, dim_j_out)                           &
     &, coeff_j_lin1   (dim_i_out, dim_j_out)                           &
     &, max_mono                                                        &
     &, min_mono

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
! Section 1.  Perform cubic Lagrange Interpolation in j direction and
!             then in i direction, for central lines, linear for outside
!             lines, and ensure all pieces are monotone.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ECMWF_MONO_QUASI_CUBIC',zhook_in,zhook_handle)
      one_sixth = 1./6.
! begin loop over levels which ends at end of subroutine

      Do k = 1, dim_k_out

        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
            phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
            phi_cubed = weight_phi (i,j,k) * phi_sq

            coeff_j_plus2(i,j) = one_sixth *                            &
     &                         (phi_cubed - weight_phi(i,j,k))
            coeff_j_plus(i,j)  = 0.5 * (phi_cubed - phi_sq              &
     &                           - 2.*weight_phi(i,j,k))
            coeff_j_zero(i,j)  = 0.5 * (phi_cubed - 2. * phi_sq -       &
     &                             weight_phi(i,j,k) + 2.)
            coeff_j_minus(i,j) = one_sixth * (phi_cubed - 3.*phi_sq +   &
     &                                    2. * weight_phi(i,j,k))
            coeff_j_lin1(i,j)  = 1.- weight_phi(i,j,k)
          End Do
        End Do

! Calculate interpolation weights in the i direction
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
            lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
            lambda_cubed = weight_lambda (i,j,k) * lambda_sq

            coeff_i_plus2(i,j) = one_sixth *                            &
     &                           (lambda_cubed - weight_lambda(i,j,k))
            coeff_i_plus(i,j)  = 0.5 * (lambda_cubed - lambda_sq        &
     &                             - 2.*weight_lambda(i,j,k))
            coeff_i_zero(i,j)  = 0.5 * (lambda_cubed - 2.*lambda_sq -   &
     &                             weight_lambda(i,j,k) + 2.)
            coeff_i_minus(i,j) = one_sixth * (lambda_cubed -            &
     &                                        3.*lambda_sq +            &
     &                                      2. * weight_lambda(i,j,k))
          End Do
        End Do

        Do n = 1, number_of_inputs
! interpolate in j at each index level for the four points needed for
! the i interpolation.
        Do index = 0, 1

          Do j = 1, dim_j_out
            Do i = 1, dim_i_out

              val_i_minus = weight_phi(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                    + coeff_j_lin1(i,j) *                         &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),k_out(i,j,k)+index,n)

                val_i     = coeff_j_plus2(i,j) *                        &
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

              max_mono = max (                                          &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+index,n),  &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+index,n) )
              min_mono = min (                                          &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+index,n),  &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+index,n) )
              If (val_i  >   max_mono ) val_i = max_mono
              If (val_i  <   min_mono ) val_i = min_mono

                val_i_plus = coeff_j_plus2(i,j) *                       &
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

              max_mono = max (                                          &
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+index,n),&
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+index,n) )
              min_mono = min (                                          &
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+index,n),&
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+index,n) )
              If (val_i_plus  >   max_mono ) val_i_plus = max_mono
              If (val_i_plus  <   min_mono ) val_i_plus = min_mono

              val_i_plus2 = weight_phi(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                    + coeff_j_lin1(i,j) *                         &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),k_out(i,j,k)+index,n)

! interpolate in i and store.

              col_data (i,j,index,n) = coeff_i_plus2(i,j)               &
     &                                   * val_i_plus2                  &
     &                                 - coeff_i_plus(i,j)              &
     &                                   * val_i_plus                   &
     &                                 + coeff_i_zero(i,j)              &
     &                                   * val_i                        &
     &                                 - coeff_i_minus(i,j)             &
     &                                   * val_i_minus
              max_mono = max (val_i, val_i_plus)
              min_mono = min (val_i, val_i_plus)
              If (col_data (i,j,index,n)  >   max_mono )                &
     &            col_data (i,j,index,n) = max_mono
              If (col_data (i,j,index,n)  <   min_mono )                &
     &            col_data (i,j,index,n) = min_mono
            End Do
          End Do

        End Do ! end loop over index

! calculate interpolation coeffs for linear
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out

            coeff_plus2(i,j) = (1.-weight_lambda(i,j,k)) *              &
     &                           (1.-weight_phi(i,j,k))
            coeff_plus(i,j)  = weight_lambda(i,j,k) *                   &
     &                           (1.-weight_phi(i,j,k))
            coeff_zero(i,j)  = (1.-weight_lambda(i,j,k)) *              &
     &                            weight_phi(i,j,k)
            coeff_minus(i,j) = weight_lambda(i,j,k) *                   &
     &                           weight_phi(i,j,k)
          End Do
        End Do

! interpolate to horizontal point doing i and j simultaneously.
! interpolate using bi-linear.
        Do index = -1, 2, 3
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out

              col_data (i,j,index,n) = coeff_plus2(i,j) *               &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+index,n)    &
     &                             + coeff_plus(i,j) *                  &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+index,n)  &
     &                             + coeff_zero(i,j) *                  &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+index,n)  &
     &                             + coeff_minus(i,j) *                 &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+index,n)


            End Do
          End Do
        End Do
        End Do ! end loop over number of inputs

! ----------------------------------------------------------------------
! Section 2.   Perform cubic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Interpolate data
        Do n = 1, number_of_inputs
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
            Data_out (i,j,k,n) = coeff_z(i,j,k,-1) *                    &
     &                         col_data (i,j,-1,n) +                    &
     &                         coeff_z(i,j,k,0) *                       &
     &                         col_data (i,j,0,n) +                     &
     &                         coeff_z(i,j,k,1) *                       &
     &                         col_data (i,j,1,n) +                     &
     &                         coeff_z(i,j,k,2) *                       &
     &                         col_data (i,j,2,n)

              max_mono = max (col_data (i,j,0,n), col_data (i,j,1,n))
              min_mono = min (col_data (i,j,0,n), col_data (i,j,1,n))
              If (Data_out (i,j,k,n)  >   max_mono )                    &
     &            Data_out (i,j,k,n) = max_mono
              If (Data_out (i,j,k,n)  <   min_mono )                    &
     &            Data_out (i,j,k,n) = min_mono
          End Do
        End Do
      End Do

! End loop over levels.
      End Do

! End of routine.
      IF (lhook) CALL dr_hook('ECMWF_MONO_QUASI_CUBIC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ECMWF_mono_quasi_cubic
