! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE BOX_SUM
!
!    NOT SUITABLE FOR SINGLE COLUMN USE
!
!
!    SUITABLE FOR ROTATED GRIDS
!
!    SYSTEM TASK:  S1 (part,extension for area mean interpolation)
!
!    PURPOSE:
!    Routine sums contributions from gridboxes for source data on a
!    regular lat-long grid to form means for gridboxes of a regular
!    lat-long grid specified as target.
!    Both grids are defined with the same pole and orientation;
!    the original data must be interpolated onto a rotated
!    grid ,if the target grid is a rotated grid, BEFORE calling this
!    routine.
!    The algorithms are general and will cope with either a finer
!    or coarser resolution source grid.
!
!    DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1
!
!
!     -------------------------------------------------------------

!    ARGUMENTS:---------------------------------------------------

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
SUBROUTINE box_sum                                                             &
  (source_row_length,source_rows,row_length,rows,                              &
  long_l,colat_t,i_l,j_t,global,boxsum,source)
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER                                                                        &
  source_row_length                                                            &
                                !IN    Number of points per row (source data)
                                !      on rotated grid if necessary
  , source_rows                                                                &
                                !IN    Number of rows of source data
                                !      on rotated grid if necessary
  , row_length                                                                 &
                                !IN    Number of points per row target area
  , rows               !IN    Number of rows of target area

INTEGER                                                                        &
  i_l(row_length+1)                                                            &
                                !IN Index of first source gridbox to overlap
                                !   with left hand side of target gridbox
  ,j_t(rows+1)      !IN Index of first source gridbox to overlap
!   bottom of target gridbox

!  N.B.I_L(I) is the first source gridbox to overlap LH side of target
!             box  I of a row
!      I_L(I+1) is the last source gridbox to overlap RH side of target
!             box  I of a row
!      J_T(J) is the first source gridbox to overlap bottom of target
!             box on row J
!      J_T(J+1) is the last source gridbox to overlap top of target
!             box on row J
!
!  REAL value of:-
!      I_L(I) is also used to measure the 'longitude' of the RHS of the
!             source gridbox
!      J_T(J) is also used to measure the 'colatitude' of the top
!             of the source gridbox


REAL                                                                           &
  source(source_row_length,source_rows)!IN  source data
REAL boxsum(row_length,rows)                                                   &
                                !OUT Sum of data on target grid
  ,    long_l(row_length +1)                                                   &
                                !IN Left longitude of gridbox (in
                                ! units of souce gridbox EW length)
  ,    colat_t(rows +1)         !IN Colatitude of bottom of gridbox (in
! units of source gridbox NS length)

LOGICAL global       !IN    true if global area required
!*---------------------------------------------------------------------

!    WORKSPACE USAGE:-------------------------------------------------
!   Define local workspace arrays:
!    1 real of length row_length
REAL                                                                           &
  ew_sum(row_length)                                                           &
                                ! summed WE source data
  ,    ew_weight(row_length)                                                   &
                                ! summed WE weights for source data
  ,    box_weight(row_length)  ! summed weights for target boxes

REAL long1, long2
INTEGER :: it1

!*---------------------------------------------------------------------

!   EXTERNAL SUBROUTINES CALLED---------------------------------------
! None
!*------------------------------------------------------------------

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!    DEFINE LOCAL VARIABLES

INTEGER i,j,i1,i2,it,j1,j2,jt,k ! loop counters

REAL weight

!  *********************************************************************
!  1.0 Sum source boxes (whole and partial) contributions to target box
!  *********************************************************************

DO j=1,rows
  j1 = j_t(j)
  j2 = j_t(j+1)

  DO i=1,row_length
    box_weight(i)=0.0
    boxsum(i,j)=0.0
  END DO

  DO jt=j1,j2

    !  *********************************************************************
    !  1.1 Sum  EW (whole and partial) contributions to target grid boxes
    !  *********************************************************************

    DO i=1,row_length
      ew_sum(i)=0.0
      ew_weight(i)=0.0
    END DO

    DO i=1,row_length
      ! Index to source point to the left of the target box (and the
      ! longitude of the right side of the source box.
      i1 = i_l(i)
      ! Index of the left coordinate for final source box (and the
      ! longitude of the right side of the source box).
      i2 = i1+MODULO(i_l(i+1)-i1,source_row_length)
      ! Longitude of the left side of target box in terms of
      ! source box grid.
      long1 = long_l(i)
      ! Longitude of the right side (left side of next target box) of
      ! target box in terms of the source box grid.
      long2 = long1 + MODULO(long_l(i+1)-long1,REAL(source_row_length))

      DO it = i1,i2
        ! The index for the source data array (not the grid).
        it1 = MODULO(it-1,source_row_length)+1
        IF ( it == i1 .AND. i2-i1 > 0) THEN
          ! Left hand side
          ! If we are the first loop and have more than one source box to
          ! use.
          weight = MODULO(REAL(it)-long1,REAL(source_row_length))
        ELSE IF (it == i2 .AND. i2-i1 > 0) THEN
          ! Right hand side
          ! If we are the last loop and have more than one source box to
          ! use.
          weight = MODULO(1.0-(REAL(it)-long2),REAL(source_row_length))
        ELSE
          ! Whole contribution
          weight = MIN(long2-long1,1.0)
        END IF
        IF (NINT(source(it1,jt)) /= NINT(rmdi)) THEN
          ew_weight(i) = ew_weight(i) + weight
          ew_sum(i)    = ew_sum(i)    + weight*source(it1,jt)
        END IF

      END DO
    END DO

    !  *********************************************************************
    !  1.3 Add summed EW  box contributions into rows J  target grid boxes
    !  *********************************************************************
    IF (jt == j1 .AND. j2-j1 > 0) THEN
      ! If we are the first source row contribution and have more than than 1
      ! source row to sum up this must be transitioning between source boxes.

      weight = (REAL(jt) - colat_t(j))
    ELSE IF (jt == j2 .AND. j2-j1 > 0) THEN
      ! If we are the last source row contribution and have more than than 1 
      ! source row to sum up this must be transitioning between source boxes.

      weight = (1.0-(REAL(jt)-colat_t(j+1)))
    ELSE
      ! Whole contributions to row J (including when only one source row is 
      ! being used for target (hence check for j2-j1 > 0 above).  We either
      ! have a whole row contribution or a proportion according to colats.

      weight = MIN(colat_t(j+1)-colat_t(j),1.0)
    END IF
    DO i=1,row_length
      boxsum(i,j)   = boxsum(i,j)   + weight*ew_sum(i)
      box_weight(i) = box_weight(i) + weight*ew_weight(i)
    END DO

  END DO

  DO i=1,row_length
    IF (box_weight(i) /= 0.0) THEN
      boxsum(i,j) = boxsum(i,j) / box_weight(i)
    ELSE
      boxsum(i,j) = rmdi
    END IF
  END DO

END DO

RETURN
END SUBROUTINE box_sum
