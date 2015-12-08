! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! SUBROUTINE BOX_BND

! Routine sets up arrays of indexes and boundaries of target grid
! boxes relative to a source grid for use  by BOX_SUM so that
! area weighted means can be calculated.

! NOT SUITABLE FOR SINGLE COLUMN USE

! SUITABLE FOR ROTATED GRIDS

! SYSTEM TASK:  S1 (part,extension for area mean interpolation)

! PURPOSE: To set up for grid-boxes on a target grid the longitude,
!          colatitude,and indexes of the overlapping source grid-
!          boxes at left hand side and top of target grid-boxes.
!          Both grids are regular lat-long with the same pole
!          and orientation. Either may be a 'p-grid' (ie with
!          half-size boxes at the poles) or a 'u-grid' (ie with
!          regular size boxes everywhere.)

!          NB The units used are "source grid-box lengths"
!          The area of the target grid_boxes in squared "source
!          grid-box units" is also returned.

! DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1

! ---------------------------------------------------------------------


! ARGUMENTS:-----------------------------------------------------------

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
SUBROUTINE box_bnd                                                             &
  (i_l,long_l,j_t,colat_t,area_box,                                            &
  row_length,rows,row_length_srce,rows_srce,                                   &
  delta_long,delta_lat,start_long,start_lat,                                   &
  delta_long_srce,delta_lat_srce,start_long_srce,start_lat_srce,               &
  igrid,igrid_srce,global)


USE ereport_mod, ONLY: ereport
USE printstatus_mod
IMPLICIT NONE

INTEGER                                                                        &
  row_length                                                                   &
                                !IN    Number of points per row target area
  , rows                                                                       &
                                !IN    Number of rows of target area
  , row_length_srce                                                            &
                                !IN    Number of points per row source area
  , rows_srce          !IN    Number of rows of source area

INTEGER                                                                        &
  i_l(row_length+1)                                                            &
                                !OUT Index of source box overlapping lhs of
                                ! target grid-box
  ,j_t(rows+1)      !OUT Index of source box overlapping bottom of
! target grid-box

REAL                                                                           &
  long_l(row_length +1)                                                        &
                                !OUT Left longitude of target grid-box (in
                                ! units of DELTA_LONG_SRCE)
  ,colat_t(rows+1)                                                             &
                                !OUT Colatitude of bottom of target grid-box
                                ! (in units of DELTA_LAT_SRCE)
  ,  area_box !OUT area of grid box in sq units of source grid

REAL                                                                           &
  delta_long                                                                   &
                                !IN   Longitude increment of target grid (deg)
  ,delta_lat                                                                   &
                                !IN   Latitude increment of target grid (deg)
  ,start_long                                                                  &
                                !IN   start longitude of centre of first grid-
                                !     box in target area
  ,start_lat                                                                   &
                                !IN   start latitude of centre of first grid-
                                !     box in target area
  ,delta_lat_srce                                                              &
                                !IN   Latitude increment of source grid
  ,delta_long_srce                                                             &
                                !IN   Longitude increment of source grid
  ,start_long_srce                                                             &
                                !IN   start longitude of centre of first grid-
                                !     box in source area
  ,start_lat_srce  !IN   start latitude of centre of first grid
!     box in source area

INTEGER                                                                        &
  igrid                                                                        &
                                !IN   Grid indicator 1=at pole,2=not at pole
  ,igrid_srce      !IN   Grid indicator 1=at pole,2=not at pole

LOGICAL global   !IN    true if global area required

!    DEFINE LOCAL VARIABLES
REAL ew_box                                                                    &
                                ! length of grid box in units of source
                                ! (DELTA_LONG_SRCE)
  ,    ns_box                                                                  &
                                ! height of grid box in units of source
                                ! (DELTA_LAT_SRCE)
  ,start_long_box                                                              &
                                ! start longitude of first grid box left edge
  ,start_colat_box                                                             &
                                ! start colatitude of first grid box top edge
  ,start_long_box_srce                                                         &
                                ! start longitude of first grid box left
                                ! edge ( source)
  ,start_colat_box_srce                                                        &
                                ! start colatitude of first grid box top
                                ! edge ( source)
  ,long_offset                                                                 &
                                ! start longitude difference
  ,colat_offset        ! start colatitude difference

INTEGER i,j ! loop counters
INTEGER :: errorstatus ! Error code for ereport
REAL p1,p2
LOGICAL lner
CHARACTER(len=*), PARAMETER  :: routinename='boxbnd'
CHARACTER(len=80)            :: cmessage
lner(p1,p2) = ((ABS(p1-p2))  >   (1.e-5*ABS(p1+p2)))

!  *********************************************************************
!  1.0 Set target gridbox length,height and area in units of source grid
!      Set start 'longitude' and 'colatitude' of first grid box
!      also in source units.
!  *********************************************************************

ew_box= delta_long/delta_long_srce
ns_box= delta_lat/delta_lat_srce
area_box= ew_box*ns_box

!  *********************************************************************
!  1.1 Set start colatitude of top and start longitude of LHS of first
!      boxes on both target and source grids.
!  *********************************************************************
IF ( printstatus >= prstatus_diag ) THEN
  WRITE(6,*) delta_lat_srce,delta_lat
  WRITE(6,*) delta_long_srce,delta_long
  WRITE(6,*) start_lat_srce,start_lat
  WRITE(6,*) start_long_srce,start_long
END IF
! The definition of colat is that it is 0.0 at the south pole and
! increasing to source_rows at the north pole.
! This is due to changing grid from north to sout
start_long_box = start_long - 0.5*delta_long
start_colat_box = (start_lat + 90.) - 0.5*delta_lat
start_long_box_srce = start_long_srce - 0.5*delta_long_srce
start_colat_box_srce = (start_lat_srce + 90.0) -0.5*delta_lat_srce

IF (global) THEN
  IF (igrid_srce == 1.AND.lner(start_lat_srce,-90.)) THEN
    cmessage = ' BOX_BND: source grid not global'
    errorstatus = 1
    CALL ereport(routinename, errorstatus, cmessage)
  END IF
  IF (igrid == 1.AND.lner(start_lat,-90.)) THEN
    cmessage = ' BOX_BND: target grid not global'
    errorstatus = 1
    CALL ereport(routinename, errorstatus, cmessage)
  END IF

ELSE
  WRITE (6,*) 'BOX_BND'
  WRITE (6,*) 'start_colat_box      ',START_COLAT_BOX
  WRITE (6,*) 'start_colat_box_srce ',START_COLAT_BOX_SRCE
  IF (lner(start_colat_box,start_colat_box_srce)) THEN
    cmessage =' BOX_BND: target area larger than source area'
    errorstatus = 1
    CALL ereport( routinename, errorstatus, cmessage)
  END IF
END IF

long_offset  =(start_long_box-start_long_box_srce)/                            &
  delta_long_srce
colat_offset = (start_colat_box - start_colat_box_srce)/                       &
  delta_lat_srce

IF (.NOT.global) THEN
  IF (long_offset <  0.0) THEN
    WRITE (6,*) ' LONG_OFFSET = ',LONG_OFFSET,' ; Reset to 0.0'
    long_offset = 0.0
  END IF
  IF (colat_offset <  0.0) THEN
    WRITE (6,*) ' COLAT_OFFSET = ',COLAT_OFFSET,' ; Reset to 0.0'
    colat_offset = 0.0
  END IF
END IF
!  *********************************************************************
!  2.0 Set grid box left longitudes, top colatitudes and indices
!  *********************************************************************

DO i=1,row_length + 1
  long_l(i) = long_offset + (i-1)*ew_box
  IF (global.AND.long_l(i) <  0.0) THEN
    long_l(i)=long_l(i)+REAL(row_length_srce)
  ELSE IF (global.AND.long_l(i) >= row_length_srce) THEN
    long_l(i)=long_l(i)-REAL(row_length_srce)
  END IF
  i_l(i) = FLOOR(long_l(i) +1.0)
END DO
IF (long_l(1) <  0.0)long_l(1)=0.0

IF ( printstatus >= prstatus_diag ) THEN
  WRITE(6,*) ' I_L'
  WRITE(6,*) i_l
  WRITE(6,*) ' LONG_L'
  WRITE(6,*) long_l
END IF

colat_t(1) = colat_offset
IF (global.AND.igrid == 1) THEN
  colat_t(1) = (0.0  - start_colat_box_srce)/delta_lat_srce
END IF
IF (colat_t(1) <  0.0)colat_t(1)=0.0
j_t(1) = FLOOR(colat_t(1) + 1.)
DO j=2,rows+1
  colat_t(j) = colat_offset  + (j-1)*ns_box
  j_t(j) = FLOOR(colat_t(j) + 1.)
END DO
!    ROWS+1 ie top boundary
IF (global) THEN
  IF (igrid == 1)colat_t(rows+1) = colat_t(rows+1)-0.5*ns_box
  j_t(rows+1) = rows_srce
ELSE
  IF (j_t(rows+1) >  rows_srce) THEN
    IF (colat_t(rows+1) >  rows_srce) THEN
      cmessage = ' BOX_BND: target area larger than source area'
      errorstatus = 1
      CALL ereport(routinename,errorstatus,cmessage)
    ELSE
      j_t(rows+1) = rows_srce
    END IF
  END IF
END IF

IF ( printstatus >= prstatus_diag ) THEN
  WRITE(6,*) ' J_T'
  WRITE(6,*) j_t
  WRITE(6,*) ' COLAT_T'
  WRITE(6,*) colat_t
END IF

RETURN
END SUBROUTINE box_bnd
