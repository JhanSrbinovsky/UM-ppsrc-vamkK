! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Given header grid information returns parameters of fixed portion of grid

SUBROUTINE ascertain_fixed_grid(len_x, len_y, x, y, origin_x, origin_y,      &
                                dx, dy, nx, ny, xstart, ystart,reqdomx, reqdomy)

USE Err_Mod, ONLY:        StatusFatal, StatusOK, StatusWarning
USE ereport_mod, ONLY: ereport
USE um_types, ONLY: real32

IMPLICIT NONE

!
! Description:
!   Finds fixed-grid regions of a variable resolution grid
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.

! Subroutine arguments  

INTEGER, INTENT(IN)  :: len_x, len_y       ! Length of input arrays
INTEGER, INTENT(IN)  :: reqdomx, reqdomy   ! Required domains in x and y
REAL, INTENT(IN)     :: x(len_x), y(len_y) ! Lat/long of variable res grid

INTEGER, INTENT(OUT) :: nx, ny             ! No. of points in x/y dir fixed grid
INTEGER, INTENT(OUT) :: xstart, ystart     ! Array posns of start of fixed grid
REAL, INTENT(OUT)    :: origin_x, origin_y ! Lat/long of first fixed res point
REAL, INTENT(OUT)    :: dx, dy             ! Fixed res grid spacing



! Local variables
INTEGER :: i                              ! DO loop variable
INTEGER, PARAMETER :: max_domains = 10    ! Maximum number of domains in one dir
REAL :: prev_spacing, next_spacing        ! Spacing on either side of current pt
REAL :: tol                               ! Very small number

! Grid parameters
INTEGER :: domainx, domainy
INTEGER :: nx_domain(max_domains), ny_domain(max_domains)
INTEGER :: xstart_domain(max_domains), ystart_domain(max_domains)
REAL    :: origin_x_domain(max_domains), origin_y_domain(max_domains)
REAL    :: dx_domain(max_domains), dy_domain(max_domains)


! Error reporting variables
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ascertain_fixed_grid"
INTEGER :: ErrorStatus                 ! Program status monitor


! End of header

! Set tolerance to be small for a 32-bit value since program generating 
! grid may have been 32-bit.
tol = EPSILON(REAL(1.0,KIND=real32))

! Initialise output variables
nx       = 0
ny       = 0
origin_x = 0
origin_y = 0
dx       = 0
dy       = 0
xstart   = 0
ystart   = 0
domainx  = 0
domainy  = 0

ErrorStatus = StatusOK


! Find the start of the fixed grid in the x direction
find_x_start:  DO i = 2, len_x-1
    prev_spacing = x(i)   - x(i-1)
    next_spacing = x(i+1) - x(i)

    IF (ABS(prev_spacing - next_spacing) < tol ) THEN

      IF(domainx > 0) THEN
        IF (ABS (next_spacing - dx_domain(domainx) ) < tol ) THEN
          CYCLE find_x_start
        END IF
      END IF

      domainx = domainx + 1
      xstart_domain(domainx) = i-1
      dx_domain(domainx) = prev_spacing
      origin_x_domain(domainx) = x(i-1)
    END IF
END DO find_x_start

WRITE(6,'("Found",I2," x domains")') domainx
DO i = 1, domainx
  WRITE(6,'(I2,1x,I5,1x,ES10.4)') i,xstart_domain(i),dx_domain(i)
END DO
   
! Failed to find fixed grid in x direction
IF (domainx == 0) THEN
  ErrorStatus = StatusFatal

  CALL EReport(RoutineName, ErrorStatus,                                     &
                            "Failed to find any start of fixed longitude grid")
END IF

IF (domainx < reqdomx) THEN
  ErrorStatus = StatusFatal

  CALL EReport(RoutineName, ErrorStatus,                                     &
                            "Failed to find required longitude domain")
END IF

! Find the start of the fixed grid in the y direction
find_y_start:  DO i = 2, len_y-1
    prev_spacing = y(i)   - y(i-1)
    next_spacing = y(i+1) - y(i)

    IF (ABS(prev_spacing - next_spacing) < tol ) THEN

      IF (domainy > 0) THEN
        IF (ABS (next_spacing - dy_domain(domainy) ) < tol ) THEN
          CYCLE find_y_start
        END IF
      END IF

      domainy = domainy + 1
      ystart_domain(domainy) = i-1
      dy_domain(domainy) = prev_spacing
      origin_y_domain(domainy) = y(i-1)
    END IF
END DO find_y_start

WRITE(6,'("Found",I2," y domains")') domainy
DO i = 1, domainy
  WRITE(6,'(I2,1x,I5,1x,ES10.4)') i,ystart_domain(i),dy_domain(i)
END DO
   
! Failed to find fixed grid in y direction
IF (domainy == 0) THEN
  ErrorStatus = StatusFatal

  CALL EReport(RoutineName, ErrorStatus,                                     &
                            "Failed to find any start of fixed latitude grid")
END IF

IF (domainy < reqdomy) THEN
  ErrorStatus = StatusFatal

  CALL EReport(RoutineName, ErrorStatus,                                     &
                            "Failed to find required latitude domain")
END IF

! Select required domain
xstart   = xstart_domain(reqdomx)
ystart   = ystart_domain(reqdomy)
dx       = dx_domain(reqdomx)
dy       = dy_domain(reqdomy)
origin_x = origin_x_domain(reqdomx)
origin_y = origin_y_domain(reqdomy)


! Find limit of fixed grid in x direction
find_x_limit: DO i = xstart+1, len_x
  prev_spacing = x(i) - x(i-1)
  IF (ABS ( prev_spacing - dx) < tol ) THEN
    CYCLE find_x_limit
  ELSE
! nx is the number of grid points with the specified fixed grid spacing
! (i.e. nx = number of spaces in grid + 1)
! When triggered, i = rightmost point + 1  (as it's doing prev_spacing)
! making the rightmost point = i - 1
! so nx = i - 1 - xstart + 1  where the + 1 is because xstart is a valid point 
! too 
    nx = i - xstart
    EXIT find_x_limit
  END IF
END DO find_x_limit

IF (nx == 0) THEN
  ErrorStatus = StatusWarning

  CALL EReport(RoutineName, ErrorStatus,                                     &
                            "Failed to find end of fixed longitude grid, assuming end")
  nx = len_x - xstart + 1
END IF
  

! Find limit of fixed grid in y direction
find_y_limit: DO i = ystart+1, len_y
  prev_spacing = y(i) - y(i-1)
  IF (ABS (prev_spacing  - dy) < tol ) THEN
    CYCLE find_y_limit
  ELSE
    ny = i - ystart
    EXIT find_y_limit
  END IF
END DO find_y_limit

IF (ny == 0) THEN
  ErrorStatus = StatusWarning

  CALL EReport(RoutineName, ErrorStatus,                                     &
                            "Failed to find end of fixed latitude grid, assuming end")
  ny = len_y - ystart + 1                         
END IF

RETURN
END SUBROUTINE ascertain_fixed_grid
