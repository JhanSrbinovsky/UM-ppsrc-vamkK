! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Cubic_Lagrange_niv

      Subroutine Cubic_Lagrange_niv(                                    &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j,                         &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, depart_level,             &
     &                          row_length_in, rows_in,                 &
     &                          lambda_rm, lambda_rp, phi_rm, phi_rp,   &
     &                          recip_lambda_m, recip_lambda_0,         &
     &                          recip_lambda_p, recip_lambda_p2,        &
     &                          recip_phi_m, recip_phi_0,               &
     &                          recip_phi_p, recip_phi_p2,              &
     &                          L_regular,                              &
     &                          model_domain,                           &
     &                          at_extremity, n_procx, n_procy,         &
     &                          global_row_length, global_rows,         &
     &                          proc_col_group, proc_row_group,         &
     &                          datastart,                              &
     &                          halo_i_out, halo_j_out,                 &
     &                          Data_out)

! Purpose:
!          Performs cubic Lagrange interpolation of the input field to a
!          set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out, but using no interpolation
!          in the vertical.
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


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_regular                                                       &
                    ! false if variable resolution
     &, at_extremity(4)    ! Indicates if this processor is at north, 
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
     &, halo_i_out                                                      &
                    ! Size of data out halo in i direction.
     &, halo_j_out                                                      &
                    ! Size of data out halo in j direction.
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
     &  Ext_Data (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j,-1:dim_k_in+2)               &
                                                          !data to be
                                                          ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

      Integer                                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
             ! lies between it and it+1
     &, depart_level (dim_i_out, dim_j_out, dim_k_out)

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

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (1-halo_i_out:dim_i_out+halo_i_out,                    &
     &            1-halo_j_out:dim_j_out+halo_j_out, dim_k_out)
                  ! data interpolated to desired locations.

! Local Variables.

      Integer                                                           &
     &  i, j, k ! Loop indices

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
     &, coeff_minus                                                     &
     &, coeff_zero                                                      &
     &, coeff_plus                                                      &
     &, coeff_plus2                                                     &
      ! x1, x3, x4, x5, x6 are temporary variables
     &, x1, x3, x4, x5, x6                                              &
     &, coeff_i_minus(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_i_zero(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_i_plus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_i_plus2(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_minus(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_zero(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_j_plus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_j_plus2(dim_i_out, dim_j_out, dim_k_out)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Perform cubic Lagrange Interpolation in j direction and
!              then in i direction.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CUBIC_LAGRANGE_NIV',zhook_in,zhook_handle)
      If ( L_regular ) Then

      one_sixth = 1./6.
! begin loop over levels which ends at end of subroutine

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)         &
!$OMP& PRIVATE(phi_sq, phi_cubed, lambda_sq, lambda_cubed,              &
!$OMP&   coeff_plus2,coeff_plus,coeff_zero,coeff_minus,                 &
!$OMP&     val_i_plus2,  val_i_plus,  val_i,       val_i_minus)         &
!$OMP& SHARED(dim_k_out,ext_data, dim_i_out,dim_j_out,weight_phi,       &
!$OMP&   one_sixth,weight_lambda,i_out,j_out,depart_level,data_out)
      Do k = 1, dim_k_out

          Do j = 1, dim_j_out
            Do i = 1, dim_i_out
! interpolate in j at each index level for the four points needed for
! the i interpolation.
              phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
              phi_cubed = weight_phi (i,j,k) * phi_sq

              coeff_plus2 = one_sixth * (phi_cubed - weight_phi(i,j,k))
              coeff_plus  = 0.5 * (phi_cubed - phi_sq                   &
     &                             - 2.*weight_phi(i,j,k))
              coeff_zero  = 0.5 * (phi_cubed - 2. * phi_sq -            &
     &                             weight_phi(i,j,k) + 2.)
              coeff_minus = one_sixth * (phi_cubed - 3. * phi_sq +      &
     &                                      2. * weight_phi(i,j,k))

              val_i_minus = coeff_plus2 *                               &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,depart_level(i,j,k))

              val_i       = coeff_plus2 *                               &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,depart_level(i,j,k))   &
     &                         - coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k))   &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k))     &
     &                         - coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,depart_level(i,j,k))

              val_i_plus  = coeff_plus2 *                               &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,depart_level(i,j,k))

              val_i_plus2 = coeff_plus2 *                               &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,depart_level(i,j,k))

! interpolate in i and store.
              lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
              lambda_cubed = weight_lambda (i,j,k) * lambda_sq

              data_out (i,j,k) = one_sixth * (lambda_cubed -            &
     &                                            weight_lambda(i,j,k)) &
     &                                         * val_i_plus2            &
     &                               - 0.5     * (lambda_cubed -        &
     &                                            lambda_sq -           &
     &                                         2.*weight_lambda(i,j,k)) &
     &                                         * val_i_plus             &
     &                               + 0.5     * (lambda_cubed -        &
     &                                            2. * lambda_sq -      &
     &                                            weight_lambda(i,j,k)  &
     &                                            + 2.)                 &
     &                                         * val_i                  &
     &                             - one_sixth * (lambda_cubed - 3. *   &
     &                                            lambda_sq + 2. *      &
     &                                            weight_lambda(i,j,k)) &
     &                                         * val_i_minus

          End Do
        End Do

      End Do  !   k = 1, dim_k_out
!$OMP END PARALLEL DO

      else  ! variable resolution

!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP& PRIVATE(i,j,k,val_i_plus2,val_i_plus,val_i,val_i_minus,x1,       &
!$OMP&   x3,x4,x5,x6)                                                   &
!$OMP& SHARED(dim_k_out,dim_i_out,dim_j_out,i_out,j_out,                &
!$OMP&   coeff_j_plus2,coeff_j_plus,coeff_j_zero,coeff_j_minus,         &
!$OMP&   coeff_i_plus2,coeff_i_plus,coeff_i_zero,coeff_i_minus,         &
!$OMP&   recip_phi_p2,recip_phi_p,recip_phi_0,recip_phi_m,              &
!$OMP&   recip_lambda_p2,recip_lambda_p,recip_lambda_0,recip_lambda_m,  &
!$OMP&   weight_phi,phi_rm,phi_rp,weight_lambda,lambda_rm,lambda_rp,    &
!$OMP&   ext_data,depart_level,data_out)
        k=1
        j=1
!$OMP DO SCHEDULE(STATIC)

        Do k = 1, dim_k_out
        Do j = 1, dim_j_out
        Do i = 1, dim_i_out




          x1 = weight_phi(i,j,k) + phi_rm(i_out(i,j,k), j_out(i,j,k))
          x3 = weight_phi(i,j,k) - 1.0
          x4 = x3 - phi_rp(i_out(i,j,k), j_out(i,j,k))
          x5 = x3 * x4
          x6 = x1 * weight_phi(i,j,k)

          coeff_j_plus2(i,j,k) = x6 * x3 *                              &
     &                           recip_phi_p2(i_out(i,j,k),j_out(i,j,k))
          coeff_j_plus(i,j,k)  = x6 * x4 *                              &
     &                           recip_phi_p(i_out(i,j,k),j_out(i,j,k))
          coeff_j_zero(i,j,k)  = x1 * x5 *                              &
     &                           recip_phi_0(i_out(i,j,k),j_out(i,j,k))
          coeff_j_minus(i,j,k) = weight_phi(i,j,k) * x5 *               &
     &                           recip_phi_m(i_out(i,j,k),j_out(i,j,k))

          x1 = weight_lambda(i,j,k) + lambda_rm(i_out(i,j,k))
          x3 = weight_lambda(i,j,k) - 1.0
          x4 = x3 - lambda_rp(i_out(i,j,k))
          x5 = x3 * x4
          x6 = x1 * weight_lambda(i,j,k)

          coeff_i_plus2(i,j,k) = x6 * x3 * recip_lambda_p2(i_out(i,j,k))
          coeff_i_plus(i,j,k)  = x6 * x4 * recip_lambda_p(i_out(i,j,k))
          coeff_i_zero(i,j,k)  = x1 * x5 * recip_lambda_0(i_out(i,j,k))
          coeff_i_minus(i,j,k) = weight_lambda(i,j,k) * x5 *            &
     &                                     recip_lambda_m(i_out(i,j,k))

        End Do  !  i loop

        END DO  !  j loop
        END DO  !  k loop

!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        Do k = 1, dim_k_out
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out

               val_i_minus       =    coeff_j_plus2(i,j,k) *            &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,depart_level(i,j,k))

               val_i             =    coeff_j_plus2(i,j,k) *            &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,depart_level(i,j,k))   &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k))   &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k))     &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,depart_level(i,j,k))

              val_i_plus         =    coeff_j_plus2(i,j,k) *            &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,depart_level(i,j,k))

              val_i_plus2        =     coeff_j_plus2(i,j,k) *           &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,depart_level(i,j,k))


           data_out (i,j,k) =  coeff_i_plus2(i,j,k) * val_i_plus2 -     &
     &                          coeff_i_plus(i,j,k) * val_i_plus +      &
     &                          coeff_i_zero(i,j,k) * val_i -           &
     &                         coeff_i_minus(i,j,k) * val_i_minus

            End Do  !  i = 1, dim_i_out
          End Do  !  j = 1, dim_j_out
        End Do  ! k = 1, dim_k_out
!$OMP END DO
!$OMP END PARALLEL

      end If ! L_regular

      IF (lhook) CALL dr_hook('CUBIC_LAGRANGE_NIV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Cubic_Lagrange_niv

