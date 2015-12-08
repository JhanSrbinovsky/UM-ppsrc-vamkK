! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_Cubic_Lagrange
!
SUBROUTINE eg_Cubic_Lagrange(                                            &
                          fld,                                           &
                          dim_i_in, dim_j_in, dim_k_in,                  &
                          dim_i_out, dim_j_out, dim_k_out,               &
                          halo_i, halo_j, number_of_inputs,              &
                          weight_lambda, weight_phi,                     &
                          s_xi1, s_xi2, t_xi1, t_xi2, q_xi1, q_xi2,      &
                          i_out, j_out, k_out,                           &
                          coeff_z, lev_ext,                              &
                          data_out)

! Purpose:
!          Performs cubic Lagrange interpolation of the input field to a
!          set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out.
!
! Method:
!          ENDGame formulation version 3.01
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

USE um_types, ONLY: integer32
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

! Array dimensions
INTEGER, INTENT(IN) :: dim_i_in     ! X dimension of input data
INTEGER, INTENT(IN) :: dim_j_in     ! Y dimension of input data
INTEGER, INTENT(IN) :: dim_k_in     ! Z dimension of input data
INTEGER, INTENT(IN) :: dim_i_out    ! X dimension of output data
INTEGER, INTENT(IN) :: dim_j_out    ! Y dimension of output data
INTEGER, INTENT(IN) :: dim_k_out    ! Z dimension of output data
INTEGER, INTENT(IN) :: halo_i       ! Halo size X
INTEGER, INTENT(IN) :: halo_j       ! Halo size Y
INTEGER, INTENT(IN) :: number_of_inputs    ! last dimension for input/output

INTEGER, INTENT(IN) :: lev_ext      ! vertical levels to loop over in coeff_z

! Input field
REAL, INTENT(IN) :: fld(1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j, &
                        -1:dim_k_in+2, number_of_inputs)

! Interpolation weights
REAL, INTENT(IN) :: weight_lambda (dim_i_out, dim_j_out, dim_k_out) 
REAL, INTENT(IN) :: weight_phi (dim_i_out, dim_j_out, dim_k_out)

REAL, INTENT(IN) :: s_xi1(2-halo_i:dim_i_in+halo_i-2)
REAL, INTENT(IN) :: s_xi2(2-halo_j:dim_j_in+halo_j-2)
REAL, INTENT(IN) :: t_xi1(2-halo_i:dim_i_in+halo_i-2)
REAL, INTENT(IN) :: t_xi2(2-halo_j:dim_j_in+halo_j-2)
REAL, INTENT(IN) :: q_xi1(-1:2,2-halo_i:dim_i_in+halo_i-2)
REAL, INTENT(IN) :: q_xi2(-1:2,2-halo_j:dim_j_in+halo_j-2)
REAL, INTENT(IN) :: coeff_z(-2:3,dim_i_out, dim_j_out, dim_k_out)

INTEGER(KIND=integer32), INTENT(IN) :: i_out (dim_i_out, dim_j_out, dim_k_out)
INTEGER(KIND=integer32), INTENT(IN) :: j_out (dim_i_out, dim_j_out, dim_k_out)
INTEGER(KIND=integer32), INTENT(IN) :: k_out (dim_i_out, dim_j_out, dim_k_out)

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT) :: data_out(dim_i_out, dim_j_out, dim_k_out, number_of_inputs)

! Local Variables.
! Loopers and indexers
INTEGER  :: i,j,k, ko,kk,ll, n, im1,i0,ip1,ip2, jm1,j0,jp1,jp2

REAL :: x, s, t      ! temporaries

REAL :: val_i_minus
REAL :: val_i
REAL :: val_i_plus
REAL :: val_i_plus2
REAL :: coeff_i_minus
REAL :: coeff_i_zero
REAL :: coeff_i_plus
REAL :: coeff_i_plus2
REAL :: coeff_j_minus
REAL :: coeff_j_zero
REAL :: coeff_j_plus
REAL :: coeff_j_plus2

! Horizontally interpolated data ready for vertical interpolation.
REAL col_data (-2:3, dim_i_out, dim_j_out )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



! ----------------------------------------------------------------------
! Section 0.   Calculate Horizontal interpolation weights
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('eg_CUBIC_LAGRANGE',zhook_in,zhook_handle)


! Interpolate data to get the data at the point for the k interpolation
!! rest of benchmark changes change results
!$OMP PARALLEL DEFAULT(SHARED)                                                 &
!$OMP& PRIVATE(i,j,k,n,val_i_plus2,val_i_plus,val_i,val_i_minus,               &
!$OMP& col_data, im1, i0, ip1, ip2, jm1, j0, jp1, jp2, ko, kk, ll,             &
!$OMP& coeff_j_plus2,coeff_j_plus,coeff_j_zero,coeff_j_minus,                  &
!$OMP& coeff_i_plus2,coeff_i_plus,coeff_i_zero,coeff_i_minus, x,s,t)
DO n = 1, number_of_inputs
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, dim_k_out

    DO j = 1, dim_j_out
      DO i = 1, dim_i_out

        i0  = i_out(i,j,k)
        im1 = i0-1
        ip1 = i0+1
        ip2 = i0+2

        j0  = j_out(i,j,k)
        jm1 = j0-1
        jp1 = j0+1
        jp2 = j0+2

        x             = weight_phi(i,j,k)
        s             = s_xi2(j0)
        t             = t_xi2(j0)

        coeff_j_minus = x*(x-1.0)*(x-t)    *q_xi2(-1,j0)
        coeff_j_zero  = (x-s)*(x-1.0)*(x-t)*q_xi2(0,j0)
        coeff_j_plus  = (x-s)*x*(x-t)      *q_xi2(1,j0)
        coeff_j_plus2 = (x-s)*x*(x-1.0)    *q_xi2(2,j0)


        x             = weight_lambda(i,j,k)
        s             = s_xi1(i0)
        t             = t_xi1(i0)

        coeff_i_minus = x*(x-1.0)*(x-t)    *q_xi1(-1,i0)
        coeff_i_zero  = (x-s)*(x-1.0)*(x-t)*q_xi1(0,i0)
        coeff_i_plus  = (x-s)*x*(x-t)      *q_xi1(1,i0)
        coeff_i_plus2 = (x-s)*x*(x-1.0)    *q_xi1(2,i0)


        ko  = k_out(i,j,k)

          DO ll = -lev_ext, lev_ext+1
            kk = ko + ll

            val_i_minus = coeff_j_plus2*fld(im1,jp2,kk,n)                &
                         +coeff_j_plus*fld(im1,jp1,kk,n)                 &
                         +coeff_j_zero*fld(im1,j0,kk,n)                  &
                         +coeff_j_minus*fld(im1,jm1,kk,n)

            val_i       = coeff_j_plus2*fld(i0,jp2,kk,n)                 &
                         +coeff_j_plus*fld(i0,jp1,kk,n)                  &
                         +coeff_j_zero*fld(i0,j0,kk,n)                   &
                         +coeff_j_minus*fld(i0,jm1,kk,n)

            val_i_plus  = coeff_j_plus2*fld(ip1,jp2,kk,n)                &
                         +coeff_j_plus*fld(ip1,jp1,kk,n)                 &
                         +coeff_j_zero*fld(ip1,j0,kk,n)                  &
                         +coeff_j_minus*fld(ip1,jm1,kk,n)

            val_i_plus2 = coeff_j_plus2*fld(ip2,jp2,kk,n)                &
                         +coeff_j_plus*fld(ip2,jp1,kk,n)                 &
                         +coeff_j_zero*fld(ip2,j0,kk,n)                  &
                         +coeff_j_minus*fld(ip2,jm1,kk,n)

            col_data(ll,i,j) = coeff_i_plus2*val_i_plus2                 &
                              +coeff_i_plus*val_i_plus                   &
                              +coeff_i_zero*val_i                        &
                              +coeff_i_minus*val_i_minus
          END DO  ! ll
        END DO ! j
      END DO ! i

! ----------------------------------------------------------------------
! Section 3.   Perform required Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Vertcal interpolate data
      SELECT CASE(lev_ext)
      CASE(0)      ! linear   vertical
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out(i,j,k,n) = coeff_z(0,i,j,k)*col_data(0,i,j) +      &
                                coeff_z(1,i,j,k)*col_data(1,i,j)

          END DO
        END DO

      CASE(1)      ! cubic   vertical
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out(i,j,k,n) = coeff_z(-1,i,j,k)*col_data(-1,i,j) +    &
                                coeff_z(0,i,j,k)*col_data(0,i,j)   +    &
                                coeff_z(1,i,j,k)*col_data(1,i,j)   +    &
                                coeff_z(2,i,j,k)*col_data(2,i,j)

          END DO
        END DO 

      CASE(2)      ! quintic vertical
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            data_out(i,j,k,n) = coeff_z(-2,i,j,k)*col_data(-2,i,j) +    &
                                coeff_z(-1,i,j,k)*col_data(-1,i,j) +    &
                                coeff_z(0,i,j,k)*col_data(0,i,j)   +    &
                                coeff_z(1,i,j,k)*col_data(1,i,j)   +    &
                                coeff_z(2,i,j,k)*col_data(2,i,j)   +    &
                                coeff_z(3,i,j,k)*col_data(3,i,j)

          END DO
        END DO
    END SELECT

  END DO ! k
!OMP END DO nowait
END DO ! n
!$OMP END PARALLEL

IF (lhook) CALL dr_hook('eg_CUBIC_LAGRANGE',zhook_out,zhook_handle)

END SUBROUTINE eg_Cubic_Lagrange
