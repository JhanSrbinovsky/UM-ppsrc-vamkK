! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE P_TO_T
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids
MODULE p_to_t_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE p_to_t (                                               &
                     ! For temperature.
! IN Field dimensions and pointers.
  row_length,rows, halo_i, halo_j                                 &
, halo_i_data, halo_j_data, levels                                &
! IN Vertical coordinate levels.
, r_theta_levels, r_rho_levels                                    &
! IN field to be interpolated.
, field_in                                                        &
! OUT Interpolated field.
, field_out                                                       &
 )


! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 90
!   THIS CODE IS WRITTEN TO UMDP3 v8.2 PROGRAMMING STANDARDS.


USE science_fixes_mod, ONLY: L_p2t_weight_fix

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

INTEGER, INTENT(IN)  ::   row_length 
INTEGER, INTENT(IN)  ::   rows
INTEGER, INTENT(IN)  ::   levels
INTEGER, INTENT(IN)  ::   halo_i
INTEGER, INTENT(IN)  ::   halo_j
INTEGER, INTENT(IN)  ::   halo_i_data
INTEGER, INTENT(IN)  ::   halo_j_data


REAL, INTENT(IN)  :: r_theta_levels (1-halo_i:row_length+halo_i,  &
                           1-halo_j:rows+halo_j,      0:levels+1)
REAL, INTENT(IN)  :: r_rho_levels (1-halo_i:row_length+halo_i,    &
                                  1-halo_j:rows+halo_j, levels+1)
REAL, INTENT(IN)  :: field_in (1-halo_i_data:row_length+halo_i_data, &
                        1-halo_j_data:rows+halo_j_data, levels+1)

! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

REAL, INTENT(OUT)  :: field_out (1-halo_i_data:row_length+halo_i_data,&
                               1-halo_j_data:rows+halo_j_data, levels)


! LOCAL VARIABLES.

INTEGER   :: i,j,k

REAL  :: weight_1 
REAL  :: weight_2
REAL  :: weight_3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Interpolate field_in (on P grid) to T grid (field_out).
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('P_TO_T',zhook_in,zhook_handle)

IF (L_p2t_weight_fix) THEN

  DO k = 1, levels
    DO j = 1-halo_j_data, rows+halo_j_data
      DO i = 1-halo_i_data, row_length+halo_i_data
        weight_1 = r_rho_levels(i,j, k+1) -                         &
                   r_rho_levels(i,j, k)
        weight_2 = r_theta_levels(i,j, k) -                         &
                   r_rho_levels(i,j, k)
        weight_3 = r_rho_levels(i,j, k+1) -                         &
                   r_theta_levels(i,j, k)
        field_out (i,j, k) =                                        &
                 weight_2/weight_1 * field_in (i,j,k+1)             &
               + weight_3/weight_1 * field_in (i,j,k)
      END DO
    END DO
  END DO
  
ELSE

  DO k = 1, levels
    DO j = 1-halo_j_data, rows+halo_j_data
      DO i = 1-halo_i_data, row_length+halo_i_data
        weight_1 = r_rho_levels(i,j, k+1) -                         &
                   r_rho_levels(i,j, k)
        weight_2 = r_theta_levels(i,j, k) -                         &
                   r_rho_levels(i,j, k)
        weight_3 = r_rho_levels(i,j, k+1) -                         &
                   r_theta_levels(i,j, k)
        field_out (i,j, k) =                                        &
                 weight_3/weight_1 * field_in (i,j,k+1)             &
               + weight_2/weight_1 * field_in (i,j,k)
      END DO
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook('P_TO_T',zhook_out,zhook_handle)
RETURN
END SUBROUTINE p_to_t
END MODULE p_to_t_mod
