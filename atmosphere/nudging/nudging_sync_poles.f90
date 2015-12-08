! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Synchronises a variable (on the theta grid) along the polar rows

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_sync_poles(        &
  proc_row_length,                    & ! Length of (LOCAL) rows in the Model
  proc_rows,                          & ! Length of (LOCAL) rows in the Model
  model_levels,                       & ! Number of levels in the Model
  halox,                              & ! Halo (EW)
  haloy,                              & ! Halo (NS)
  variable,                           & ! Variable to synchronise
  sinangle_min,                       & ! Minimum angle in the PE
  sinangle_max,                       & ! Maximum angle in the PE
  debug)                                ! Debug flag


USE conversions_mod, ONLY: pi
USE nudging_control                     ! use standard nudging swictehs


USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE UM_ParVars
IMPLICIT NONE


INTEGER, INTENT(IN) :: proc_row_length    ! local row length
INTEGER, INTENT(IN) :: proc_rows          ! local rows
INTEGER, INTENT(IN) :: model_levels       ! model levels
INTEGER, INTENT(IN) :: halox              ! Halo (EW)
INTEGER, INTENT(IN) :: haloy              ! Halo (NS)

REAL, INTENT(INOUT)  :: variable(                                  &
 (1-halox):(proc_row_length+halox),                                &
 (1-haloy):(proc_rows+haloy),                                      &
  1:model_levels )                        ! Model variable

REAL, INTENT(IN)   :: sinangle_min      ! Minimum angle in the PE
REAL, INTENT(IN)   :: sinangle_max      ! Maximum angle in the PE
INTEGER, INTENT(IN):: debug             ! debug flag

REAL   :: arctic_sum                                              &
 (proc_row_length, model_levels)       ! sum arctic values
REAL   :: antarctic_sum                                           &
 (proc_row_length, model_levels)       ! sum antarctic values
REAL   :: arctic_ave                                              &
 (model_levels)                        ! average arctic values
REAL   :: antarctic_ave                                           &
 (model_levels)                        ! average antarctic values

INTEGER   :: ii, jj                    ! Loop vraiables
INTEGER   :: dummy                     ! Dummy argument for gsum
INTEGER   :: polar_latitude            ! Latitude of Pole
LOGICAL   :: polar                     ! Are we at the Pole
LOGICAL   :: arctic                    ! Artic (T)/Antarctic (F)

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_SYNC_POLES',zhook_in,zhook_handle)

!*******************************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*) 'PE', mype,                                        &
  ': NUDGING_SYNC_POLES: Entering Routine'
END IF

! If the sine of the latitude is equal to +/-1 then we are at the poles
! If this is the case force synchronisation for the latitude specified
!***********************************************
IF(sinangle_max == 1.0) THEN

! Say where we are
  IF(debug > 10) THEN
    WRITE(OUT,*) ' NUDGING_SYNC_POLES:  Maximum angle is at Pi/2'
  END IF

! Set flags that let us knwo we are in the arctic
  polar=.TRUE.
  arctic=.TRUE.

! Our pole has occured at the end of our loop
  polar_latitude = proc_rows

! Sum the values to the arctic sum and not the antarctic one
  CALL nudging_sync_pole_sum(arctic, polar_latitude)
  antarctic_sum = 0

!***********************************************
ELSE IF(sinangle_max == -1.0) THEN

! Say where we are
  IF(debug > 10) THEN
    WRITE(OUT,*) 'NUDGING_SYNC_POLES:  Maximum angle is at -Pi/2'
  END IF

! Set flags that let us knwo we are in the antarctic
  polar=.TRUE.
  arctic=.FALSE.

! Our pole has occured at the end of our loop
  polar_latitude = proc_rows

! Sum the values to the antarctic sum and not the arctic one
 CALL nudging_sync_pole_sum(arctic, polar_latitude)
 arctic_sum = 0

!***********************************************
ELSE IF(sinangle_min == 1.0) THEN

! Say where we are
  IF(debug > 10) THEN
    WRITE(OUT,*) 'NUDGING_SYNC_POLES:  Minimum angle is at Pi/2'
  END IF

! Set flags that let us knwo we are in the arctic
  polar=.TRUE.
  arctic=.TRUE.

! Our pole has occured at the start of our loop
  polar_latitude = 1

! Sum the values to the arctic sum and not the antarctic one
 CALL nudging_sync_pole_sum(arctic, polar_latitude)
 antarctic_sum = 0

!***********************************************
ELSE IF(sinangle_min == -1.0)  THEN

! Say where we are
  IF(debug > 10) THEN
    WRITE(OUT,*) ' NUDGING_SYNC_POLES:  Minimum angle is at -Pi/2'
  END IF

! Set flags that let us knwo we are in the antarctic
  polar=.TRUE.
  arctic=.FALSE.

! Our pole has occured at the start of our loop
  polar_latitude = 1

! Sum the values to the antarctic sum and not the arctic one
  CALL nudging_sync_pole_sum(arctic, polar_latitude)
  arctic_sum = 0

!***********************************************
ELSE
! Say we are not at the Poles
  IF(debug > 10) THEN
    WRITE(OUT,*) ' NUDGING_SYNC_POLES:  Not at the Poles'
  END IF

! Set flags that let us know that we are away from the poles
  polar=.FALSE.

! Therefore this makes no contribution to the sum over polar values
 arctic_sum = 0.0
  antarctic_sum = 0.0

END IF      ! sin latitude

! Standardise the processors
CALL gc_ssync(nproc_x,dummy)

IF(debug > 10) THEN
  WRITE(OUT,*)' NUDGING_SYNC_POLES:  About to call GC_RSUM'
END IF

! Sum up arctic values for the purpose of averaging
CALL gc_rsum (                       &
  proc_row_length*model_levels,      &  ! Number of points run over
  nproc,                             &  ! Number of processors
  dummy,                             &  ! communication flag
  arctic_sum                         &  ! Temporary array
)

IF(debug > 10) THEN
  WRITE(OUT,*) ': NUDGING_SYNC_POLES:  About to call GC_SSYNC'
END IF

! Standardise the processors
CALL gc_ssync(nproc_x,dummy)

IF(debug > 10) THEN
  WRITE(OUT,*) ' NUDGING_SYNC_POLES:  About to call GC_RSUM (II)'
END IF

! Sum up antarctic values for the purpose of averaging
CALL gc_rsum (                        &
  proc_row_length*model_levels,       &  ! Number of points run over
  nproc,                              &  ! Number of processors
  dummy,                              &  ! communication flag
  antarctic_sum                       &  ! Temporary array
)

IF(debug > 10) THEN
  WRITE(OUT,*) 'NUDGING_SYNC_POLES:  About to call GC_SSYNC (II)'
END IF

! Standardise the processors
CALL gc_ssync(nproc_x,dummy)

! Are we in polar regions
IF(polar) THEN

! If we are in the  arctic / antarctic
  IF(arctic) THEN

    IF(debug > 10) THEN
      WRITE(OUT,*)' NUDGING_SYNC_POLES:  Greenland Here we come!'
    END IF

! Initialise average sum to zero
    arctic_ave = 0

! Loop over all the levels
    DO ii=1, model_levels

! Loop round colums summing up the processor sums
      DO jj=1, proc_row_length
        arctic_ave(ii) = arctic_ave(ii)                              &
                         + arctic_sum(jj, ii)
      END DO ! Columns

! Calculate average by dividing by no. processors*no. columns
      arctic_ave(ii) = arctic_ave(ii)                                &
                      / (nproc_x*proc_row_length)

! Loop round colums setting variable to be the average
      DO jj=1, proc_row_length
        variable(jj, polar_latitude, ii) = arctic_ave(ii)
      END DO ! Columns

    END DO ! Levels

! If in Antarctic
  ELSE

    IF(debug > 10) THEN
      WRITE(OUT,*) ' NUDGING_SYNC_POLES:  Antarctica Here we come!'
    END IF

! Initialise average sum to zero
    antarctic_ave = 0

! Loop over all the levels
    DO ii=1, model_levels

! Loop round colums summing up the processor sums
      DO jj=1, proc_row_length
        antarctic_ave(ii) = antarctic_ave(ii)                          &
                            + antarctic_sum(jj, ii)

        IF(debug > 10.AND.ii == 3) THEN
          WRITE(OUT,*) ': NUDGING_SYNC_POLES:  Antarctica Summed=',     &
           antarctic_sum(jj,ii)
        END IF

      END DO ! Columns

! Calculate average by dividing by no. processors*no. columns
      antarctic_ave(ii) = antarctic_ave(ii)                            &
                         / (nproc_x*proc_row_length)

! Loop round colums setting variable to be the average
      DO jj=1, proc_row_length
        variable(jj, polar_latitude, ii) = antarctic_ave(ii)
      END DO ! Columns

    END DO ! Levels

  END IF ! N Pole/ S Pole

! If not at the poles
ELSE

  IF(debug > 10) THEN
    WRITE(OUT,*) ' NUDGING_SYNC_POLES:  Not at the Poles, so do nothing'
  END IF

END IF ! Polar

IF(debug > 10) THEN
  WRITE(OUT,*)'NUDGING_SYNC_POLES: Sums equal',                        &
   arctic_sum(1,13), antarctic_sum(1,13),                              &
   arctic_ave(13), antarctic_ave(13),                                  &
   variable(1,1,13)
END IF

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*) ' NUDGING_SYNC_POLES: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_SYNC_POLES',zhook_out,zhook_handle)

RETURN

CONTAINS
! #############################################################################
! Takes latitude point as input and then loops round levels
! At each level adds up all the values of the variable at the poles

SUBROUTINE nudging_sync_pole_sum(arctic,n)

IMPLICIT NONE
INTEGER, INTENT(IN) :: n  ! Choice of latitude
LOGICAL, INTENT(IN) :: arctic ! N Pole/ S Pole
INTEGER             :: i, j  ! Loop variables
REAL(KIND=jprb)     :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_SYNC_POLES_SUM',zhook_in,zhook_handle)

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*) ' NUDGING_SYNC_POLE_SUM: Entering Routine'
END IF

DO i=1, model_levels
  DO j=1, proc_row_length

   IF(arctic) THEN
     arctic_sum(j, i) = variable(j, n, i)
   ELSE
     antarctic_sum(j, i) = variable(j, n, i)

     IF(debug > 10.AND.i == 3) THEN
       WRITE(OUT,*) ' NUDGING_SYNC_POLES:  Antarctica =',              &
                    antarctic_sum(j,i)
      END IF
   END IF

  END DO ! Longitude
END DO ! Levels

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*) ' NUDGING_SYNC_POLE_SUM: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_SYNC_POLES_SUM',zhook_out,zhook_handle)

RETURN

END SUBROUTINE nudging_sync_pole_sum
! ############################################################################

END SUBROUTINE nudging_sync_poles

