! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Constructs simple 2D relaxation parameter for nudging

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_NUDGING_CONTROL.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------

SUBROUTINE nudging_call_relax( &
  proc_row_length_min,         &  ! Min. column 
  proc_row_length_max,         &  ! Max. column 
  proc_rows_min,               &  ! Minimum local row
  proc_rows_max,               &  ! Maximum local row
  model_levels,                &  ! Number of levels
  varname,                     &  ! Variable name
  timestep,                    &  ! Timestep size
  tropopause_level,            &  ! Tropopause level
  relaxation_parameter,        &  ! Relaxation Parameter
  debug)                          ! Miscellanea

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!*******************************************************************

INTEGER, INTENT(IN) :: proc_row_length_min   ! Minimum local column
INTEGER, INTENT(IN) :: proc_row_length_max   ! Maximum local column
INTEGER, INTENT(IN) :: proc_rows_min         ! Minimum local row
INTEGER, INTENT(IN) :: proc_rows_max         ! Maximum local row
INTEGER, INTENT(IN) :: model_levels          ! Number of levels
CHARACTER(LEN=*)        :: varname               ! Variable name
REAL, INTENT(IN)    :: timestep              ! Timestep size

INTEGER, INTENT(IN)   :: tropopause_level(                          &
 proc_row_length_min:proc_row_length_max,                        &
 proc_rows_min:proc_rows_max)                ! Tropopause level

REAL, INTENT(OUT)   :: relaxation_parameter (                    &
 proc_row_length_min:proc_row_length_max,                        &
 proc_rows_min:proc_rows_max,1:model_levels) ! Relaxation Parameter

INTEGER, INTENT(IN) :: debug                 ! Debug flag

REAL                :: k_factor (                                &
 proc_row_length_min:proc_row_length_max,                        &
 proc_rows_min:proc_rows_max,1:model_levels)  ! k factor

INTEGER             :: i, j, k               ! Loop variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_CALL_RELAX',zhook_in,zhook_handle)

!**********************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ': NUDGING_CALL_RELAX: Entered Routine', timestep
END IF

! initialise variables
! relaxation parameter sets magnitude and k factor gives height dependence
relaxation_parameter(:,:,:) = 0.00
k_factor(:,:,:)             = 1.00

! Loop over rows and levels setting relaxation parameter values
DO i = proc_rows_min, proc_rows_max
  DO j = proc_row_length_min , proc_row_length_max
    DO k = ndg_lev_bottom, ndg_lev_top

! Calculate height dependence (k factor)
! Start at chosen top and cycle to chosen bottom
! Have linear turn ons at top and bottom
! All parameters determined by the nudging control module
      IF(k < (ndg_lev_bottom + ndg_on_lev_bottom)) THEN
        k_factor(j,i,k) = k_factor(j,i,k)                        &
        *(1./(ndg_on_lev_bottom))                                &
        *(k-ndg_lev_bottom)
      ELSE
        IF(k > (ndg_lev_top-ndg_on_lev_top)) THEN
          k_factor(j,i,k) = k_factor(j,i,k)                      &
          *(1./(ndg_on_lev_top))                                 &
          *(ndg_lev_top-k)
        ELSE
          k_factor(j,i,k) = 1.0
        END IF
      END IF

! Use different magnitudes of nudging in troposphere and stratosphere  
      IF(k > (tropopause_level(j,i))) THEN 
        k_factor(j,i,k) = k_factor(j,i,k)                          & 
        * ndg_strat_fac 
      ELSE IF(k>(tropopause_level(j,i)-ndg_on_lev_top)) THEN 
        k_factor(j,i,k) = k_factor(j,i,k)                          & 
        *(((1./(ndg_on_lev_top))                                   & 
        *(tropopause_level(j,i)-k)) + ndg_strat_fac) 
      END IF 

! This calculation is a little sloppy
! If turn on factor lies outside our bounds then cutoff
      IF(k_factor(j,i,k) > 1.0) k_factor(j,i,k) = 1.0
      IF(k_factor(j,i,k) < 0.0) k_factor(j,i,k) = 0.0

! Set to variable specific constant
      IF(varname == u_name) THEN
          relaxation_parameter(j,i,k) = ndg_relax_uvalue
      ELSE IF(varname == v_name) THEN
          relaxation_parameter(j,i,k) = ndg_relax_vvalue
      ELSE IF(varname == temp_name) THEN
          relaxation_parameter(j,i,k) = ndg_relax_tvalue
      ELSE
          relaxation_parameter(j,i,k) = 0.0
      END IF

! combine all factors to produce relaxation parameter
      relaxation_parameter(j,i,k) = relaxation_parameter(j,i,k)    &
                                 *timestep                         &
                                 *k_factor(j,i,k)
   END DO ! loop over levels
  END DO   ! loop over columns
END DO       ! loop over rows

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                     &
  ' NUDGING_CALL_RELAX: Entered Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_CALL_RELAX',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_call_relax
