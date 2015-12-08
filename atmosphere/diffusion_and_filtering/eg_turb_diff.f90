! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_turb_diff

      Subroutine eg_turb_diff(                                          &
                              field, fld_type,                          &
                              row_length, rows,                         &
                              in_rows, o_rows,                          &
                              model_levels, levels,                     &
                              off_x, off_y, halo_i, halo_j,             &
                              off_i, off_j, halo_x, halo_y,             &
                              xi1, xi1_off, xi2, xi2_off,               &
                              csxi2, csxi2_off,                         &
                              coeff_lam, coeff_phi,                     &
                              delta_z, recip_r_squared_delz,            &
                              j_start, j_stop,                          &
                              L_swap, field_star)

      USE timestep_mod, ONLY: timestep

! Purpose:
!          Turbulent diffusion of a field
!          Based on conservative diffusion operator.
!          v-at-poles version
!
! Method:
!          Is described in ;
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      LOGICAL ::  L_swap   ! swap_bounds logical 

      INTEGER                                                           &
     &  row_length                                                      &
                           ! number of point on a row.
     &, rows                                                            &
                           ! number of rows.
     &, in_rows                                                         &
                           ! number of field rows.
     &, o_rows                                                          &
                           ! number of offset rows
     &, model_levels                                                    &
                           ! number of model levels.
     &, levels                                                          &
                           ! number of levels to process
     &, off_x                                                           &
                           ! Size of small halo in i
     &, off_y                                                           &
                           ! Size of small halo in j.
     &, halo_i                                                          &
                           ! Size of halo in i direction.
     &, halo_j                                                          &
                           ! Size of halo in j direction.
     &, off_i                                                           &
                           ! Size of small halo for field_star
     &, off_j                                                           &
                           ! Size of small halo for field_star
     &, halo_x                                                          &
                           ! Size i-halo for field.
     &, halo_y                                                          &
                           ! Size j-halo for field
     &, j_start, j_stop                                                 &
                           ! j start/end points                         &
     &, fld_type      ! field type.

REAL, INTENT(IN) ::                                               &
  xi1(1-halo_i:row_length+halo_i)                                 &
, xi1_off(1-halo_i:row_length+halo_i)                             &
, xi2(1-halo_j:in_rows+halo_j)                                    &
, xi2_off(1-halo_j:o_rows+halo_j)                                 &
, csxi2(1-off_y:in_rows+off_y)                                    &
, csxi2_off(1-off_y:o_rows+off_y)

REAL, INTENT(IN) ::                                                     &
  coeff_lam(0:row_length+1, 0:rows+1, levels)                         &
, coeff_phi(0:row_length+1, 0:rows+1, levels)                         &
, delta_z(  0:row_length+1, 0:rows+1, levels )                        &
, recip_r_squared_delz(row_length, rows, levels )

      Real                                                              &
     &  field (1-halo_x:row_length+halo_x,                              &
     &         1-halo_y:in_rows+halo_y, model_levels )

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  field_star (1-off_i:row_length+off_i,                           &
     &              1-off_j:in_rows+off_j, model_levels )

! Local Variables.

      Integer                                                           &
     &  i, j, k     ! Loop indices

! Local arrays

      Real                                                              &
     &  temp(1-off_x:row_length+off_x, 1-off_y:rows)                    &
     &, diff_term(row_length, in_rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set values and calculate delta_z
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_TURB_DIFF',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 2.0  Horizontal Diffusion
! ----------------------------------------------------------------------
!
! "conservative diffusion" diffuses dzQ rather than Q.  Use 
! dz*(dQ/d_lambda) instead of 1/d_lambda * (dzQ) (and similarly for phi)
! to avoid instabilities ocurring in the presence of orography due to 
! volume differences between neighbouring gridboxes.

      Do k = 1, levels

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------
        Do j = j_start, j_stop
          Do i = 0, row_length
            temp(i,j) = ( field(i+1,j,k) - field(i,j,k) ) *             &
                                           delta_z(i,j,k) *             &
                                           coeff_lam(i,j,k) /           &
                                   ( csxi2_off(j) * csxi2_off(j) *      &
                                    (xi1_off(i+1) - xi1_off(i)) )    
         End Do

          Do i = 1, row_length
            diff_term(i,j) = ( temp(i,j) - temp(i-1,j) ) /              &
                               ( xi1(i) - xi1(i-1) )
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

        Do j = j_start-1, j_stop
          Do i = 1, row_length
            temp(i,j) = ( field(i,j+1,k)  -  field(i,j,k) ) *           &
                                             delta_z(i,j,k) *           &
                                           coeff_phi(i,j,k) *           &
                                           csxi2_off(j)
          End Do
        End Do

        Do j = j_start, j_stop
          Do i = 1, row_length
            diff_term(i,j) = diff_term(i,j) +                           &
                             ( temp(i,j) - temp(i,j-1) ) /              &
                             ( csxi2_off(j) * (xi2(j) - xi2(j-1)) )
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do j = j_start, j_stop
          Do i = 1, row_length
              field_star(i,j,k) = field_star(i,j,k) + timestep *        &
                                     recip_r_squared_delz(i,j,k) *      &
                                                diff_term(i,j)
          End Do
        End Do

      End Do ! k = 1, levels
      
! DEPENDS ON: swap_bounds
      IF (off_i > 0 )CALL swap_bounds(                                  &
                               field_star, row_length, in_rows, levels, &
                               off_i, off_j, fld_type, L_swap)

      IF (lhook) CALL dr_hook('EG_TURB_DIFF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_turb_diff
