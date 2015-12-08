! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE TRACER_MASS
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Tracer Advection
      SUBROUTINE TRACER_MASS(                                           &
     &                      halo_i, halo_j, off_x, off_y                &
     &,                     halo_i_field,halo_j_field                   &
     &,                     mype                                        &
     &,                     proc_all_group                              &
     &,                     row_length, rows, n_rows                    &
     &,                     model_levels,wet_model_levels,field_levels  &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     FV_cos_theta_latitude                       &
     &,                     delta_lambda,delta_phi                      &
     &,                     rho_r2, q, qcl,qcf,field1,field2,field3     &
     &,                     exner_theta_levels,cprint)

!  PURPOSE : Check tracer mass at a given point in code
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      USE global_2d_sums_mod, ONLY : global_2d_sums

      IMPLICIT NONE

! constants required for calculations
! Required for interpolation

!----------------------------------------------------------------------
! INPUT variables
!----------------------------------------------------------------------
! MPP & halo info

      Integer off_x, off_y, halo_i, halo_j   ! halo info
      Integer halo_i_field, halo_j_field   ! halo info for field array
      Integer                                                           &
     & mype               ! processor number

      Integer :: proc_all_group            ! group of all processors

! input model dimensional info
      Integer                                                           &
     &   row_length                                                     &
                             ! number of points per row
     &,  rows                                                           &
                             ! number of rows on p grid
     &,  n_rows                                                         &
                             ! number of rows on v grid
     &,  model_levels                                                   &
                             ! number of model levels
     &,  field_levels                                                   &
                             ! number of field levels
     &,  wet_model_levels    ! number of wet model levels

! Input radius of model levels
      REAL                                                              &
     &  R_THETA_LEVELS(1-halo_i:row_length+halo_i,                      &
     &            1-halo_j:rows+halo_j,0:model_levels)                  &
     &, R_RHO_LEVELS(1-halo_i:row_length+halo_i,                        &
     &            1-halo_j:rows+halo_j,model_levels)

! Input trigonometric functions (required by polar_wind)
      Real                                                              &
     & FV_cos_theta_latitude(1-off_x:row_length+off_x,                  &
     &                        1-off_y:rows+off_y)

! Input model grid spacing in radians

      REAL delta_lambda,delta_phi

!Input model data
      REAL                                                              &
     &  rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &        model_levels)                                             &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_model_levels)                                           &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      wet_model_levels)                                           &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      wet_model_levels)                                           &
     &, field1(1-halo_i_field:row_length+halo_i_field,                  &
     &         1-halo_j_field:rows+halo_j_field,                        &
     &         field_levels)                                            &
     &, field2(1-halo_i_field:row_length+halo_i_field,                  &
     &         1-halo_j_field:rows+halo_j_field,                        &
     &         field_levels)                                            &
     &, field3(1-halo_i_field:row_length+halo_i_field,                  &
     &         1-halo_j_field:rows+halo_j_field,                        &
     &         field_levels)                                            &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)
      CHARACTER(LEN=30) cprint

!In/OUT  Vertical integrals plus info on rho_dry

      INTEGER N_SUMS
      parameter (n_sums=5)
      REAL                                                              &
     &  vert_int_array(row_length,rows,n_sums)    ! vertical integrals

!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------
! pointers for sum_array
      Integer                                                           &
     & ip_q,ip_qcl,ip_qcf                                               &
     &,ip_dry                                                           &
     &,ip_wet

      Parameter (ip_q=1,ip_qcl=2,ip_qcf=3,                              &
     &           ip_dry=4, ip_wet=5)

      Real                                                              &
     &  weight1, weight2, weight3                                       &
                                  ! weights
     &, tempd, tempw, ww2, ww1                                          &
     &, tot_q                                                           &
     &, factor                         ! grid resolution factor

      Real                                                              &
     &  rho_dry(row_length,rows,model_levels)                           &
                                                ! rho dry x r^2
     &, delr_rho(row_length,rows,model_levels)                          &
                                                ! dr
     &, field1_rho(row_length,rows,field_levels)                        &
     &, field2_rho(row_length,rows,field_levels)                        &
     &, field3_rho(row_length,rows,field_levels)                        &
     &,SUM_RESULTS(n_sums)                  ! sum of SUM_ARRAY


      INTEGER I, J, K                ! LOOP COUNTER

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!----------------------------------------------------------------------
! zero summation arrays
!----------------------------------------------------------------------
       IF (lhook) CALL dr_hook('TRACER_MASS',zhook_in,zhook_handle)
      Do k=1,n_sums
        sum_results(k)=0.0
      Enddo
      Do k=1,n_sums
        DO j = 1, rows
          DO I=1,row_length
            vert_int_array(I,J,k) = 0.0
          enddo
        enddo
      enddo

!----------------------------------------------------------------------
! CALCULATE layer thickness for rho layers and calculate rho
!----------------------------------------------------------------------
      DO K=1,model_levels
        Do j = 1, rows
          DO I=1,row_length
            DELR_RHO(I,J,K) = R_THETA_LEVELS(I,J,K) -                   &
     &                        R_THETA_LEVELS(I,J,K-1)
          END DO
        END DO
      END DO

!----------------------------------------------------------------------
! convert rho to rho dry in the same way as done in dynamics - Flux_rho
! Note now uses linear vertical intepolation to be consistent with
! new dynamics code at UM 5.1.
!  Use start of timestep rho and q (or end of timestep values)
!----------------------------------------------------------------------
      k = 1
        Do j = 1, rows
          Do i = 1, row_length
            rho_dry(i,j,k) = rho_r2(i,j,k) *                            &
     &                   (1. - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
          End Do
        End Do

      Do k = 2, wet_model_levels
        Do j = 1, rows
          Do i = 1, row_length
! linear interpolation weights
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)

            tempd = ( weight2 *                                         &
     &                (1. - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)) +       &
     &               weight1 *                                          &
     &                (1. - q(i,j,k-1) - qcl(i,j,k-1) - qcf(i,j,k-1)) ) &
     &              / weight3

            rho_dry(i,j,k) = rho_r2(i,j,k) * tempd
          End Do
        End Do
      End Do
      k =wet_model_levels+1
      if(k <= model_levels)then
        Do j = 1, rows
          Do i = 1, row_length
! linear interpolation weights
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)

            tempd = ( weight2  +                                        &
     &                weight1 *                                         &
     &     (1. - q(i,j,k-1) - qcl(i,j,k-1) - qcf(i,j,k-1)) ) / weight3
            rho_dry(i,j,k) = rho_r2(i,j,k) * tempd
          End Do
        End Do
      endif

      Do k = wet_model_levels+2,model_levels
        Do j = 1, rows
          Do i = 1, row_length
            rho_dry(i,j,k) = rho_r2(i,j,k)
          End Do
        End Do
      End Do
!-------------------------------------------------------------------
! Intepolate T, w & q to rho points
! Note needs to remain consistent with intepolation used in dynamics.
! Using linear interpolation
!
!                      K               for rho
!      K                          K-1  for theta
!      X<--- w1------->X<-- w2--->X
!       <----------w3------------>
!-------------------------------------------------------------------
      k=1
        Do j = 1, rows
          DO I=1,row_length
! assume bottom rho level value equal to bottom theta level value
            field1_rho(i,j,k)   = field1(i,j,k)
            field2_rho(i,j,k)   = field2(i,j,k)
            field3_rho(i,j,k)   = field3(i,j,k)
          END DO
        END DO

      DO K=2,field_levels
        Do j = 1, rows
          DO I=1,row_length
            weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
            weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
            weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
            ww1 = weight1/weight3
            ww2 = weight2/weight3
            field1_rho (i,j, K)  = ww2 * field1(i,j,K)                  &
     &                               + ww1 * field1(i,j,K-1)
            field2_rho (i,j, K)  = ww2 * field2(i,j,K)                  &
     &                               + ww1 * field2(i,j,K-1)
            field3_rho (i,j, K)  = ww2 * field3(i,j,K)                  &
     &                               + ww1 * field3(i,j,K-1)
          END DO
        END DO
      END DO

! vertical integrals
!--------------------
      DO K=1,model_levels
        DO j = 1, rows
          DO I=1,row_length

            tempw = RHO_r2(I,j,k)*DELR_RHO(I,j,k)  ! wet mass
            tempd = RHO_dry(I,j,k)*DELR_RHO(I,j,k)  ! dry mass

! dry mass
            vert_int_array(I,J,ip_dry) =                                &
     &                    vert_int_array(i,j,ip_dry) + tempd
! wet  mass
            vert_int_array(I,J,ip_wet) =                                &
     &                    vert_int_array(i,j,ip_wet) + tempw

          enddo
        enddo
      enddo

! (b) fields on field levels only

      DO K=1,field_levels
        DO j = 1, rows
          DO I=1,row_length
            tempd = RHO_dry(I,j,k)*DELR_RHO(I,j,k)  ! dry mass
! q*rho
            vert_int_array(I,J,ip_q) = vert_int_array(i,j,ip_q) +       &
     &                                      field1_rho(i,j,k)*tempd
            vert_int_array(I,J,ip_qcl) = vert_int_array(i,j,ip_qcl) +   &
     &                                      field2_rho(i,j,k)*tempd
            vert_int_array(I,J,ip_qcf) = vert_int_array(i,j,ip_qcf) +   &
     &                                      field3_rho(i,j,k)*tempd

          END DO
        END DO
      END DO
!----------------------------------------------------------------------
      Do k = 1,n_sums
        do j = 1,rows
          do i = 1,row_length
             vert_int_array(i,j,k) = vert_int_array(i,j,k)*             &
     &                                FV_cos_theta_latitude(i,j)
          enddo
        enddo
      enddo

      CALL global_2d_sums(vert_int_array, row_length, rows, 0, 0,       &
                          n_sums, sum_results, proc_all_group)

      factor=delta_lambda*delta_phi   ! already in radians


      if (mype == 0) then
        write(6,*)'position of check ',cprint
        write(6,*)'field 1       ',sum_results(ip_q)*factor
        write(6,*)'field 2       ',sum_results(ip_qcl)*factor
        write(6,*)'field 3       ',sum_results(ip_qcf)*factor
        write(6,*)'dry mass       ',sum_results(ip_dry)*factor
        write(6,*)'wet mass       ',sum_results(ip_wet)*factor
      endif

      IF (lhook) CALL dr_hook('TRACER_MASS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TRACER_MASS
