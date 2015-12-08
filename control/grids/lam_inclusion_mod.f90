! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lam_inclusion_mod

! Description:
!   This checks that the target grid is within the source grid.  Used my makebc
!   and the reconfiguration to make sure interpolations can be considered safe.
!
! Code Owner: See Unified Model Code Owners HTML page
! 
! This file belongs in section: Grids
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

IMPLICIT NONE

CONTAINS
SUBROUTINE lam_inclusion(src_rows , src_row_len , src_halo_x , src_halo_y ,   &
                         src_lambda, src_phi, src_pole_lat, src_pole_lon,     &
                         targ_rows, targ_row_len, targ_halo_x, targ_halo_y,   &
                         targ_lambda, targ_phi, targ_phi_min, targ_phi_max,   &
                         targ_pole_lat, targ_pole_lon,                        &
                         l_same_rotation, l_src_rotated, grid_tol_arg)

USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE lltoeq_mod,       ONLY: lltoeq
USE eqtoll_mod,       ONLY: eqtoll
USE ereport_mod,      ONLY: ereport
USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

INTEGER, INTENT(IN) :: src_rows
INTEGER, INTENT(IN) :: src_row_len
INTEGER, INTENT(IN) :: src_halo_x
INTEGER, INTENT(IN) :: src_halo_y
REAL   , INTENT(IN) :: src_lambda(1-src_halo_x: src_row_len + src_halo_x)
REAL   , INTENT(IN) :: src_phi   (1-src_halo_y: src_rows + src_halo_y)
REAL   , INTENT(IN) :: src_pole_lat
REAL   , INTENT(IN) :: src_pole_lon

INTEGER, INTENT(IN) :: targ_rows
INTEGER, INTENT(IN) :: targ_row_len
INTEGER, INTENT(IN) :: targ_halo_x
INTEGER, INTENT(IN) :: targ_halo_y
REAL   , INTENT(IN) :: targ_lambda(1-targ_halo_x: targ_row_len + targ_halo_x)
REAL   , INTENT(IN) :: targ_phi   (1-targ_halo_y: targ_rows + targ_halo_y)
REAL   , INTENT(IN) :: targ_phi_min
REAL   , INTENT(IN) :: targ_phi_max
REAL   , INTENT(IN) :: targ_pole_lat
REAL   , INTENT(IN) :: targ_pole_lon

LOGICAL, INTENT(IN) :: l_same_rotation
LOGICAL, INTENT(IN) :: l_src_rotated

! Optional arguments
REAL, OPTIONAL, INTENT(IN) :: grid_tol_arg

! Local
INTEGER :: ipt
INTEGER :: j
INTEGER :: errorstatus
INTEGER :: num_pts

REAL    :: src_phi_min
REAL    :: src_phi_max
REAL    :: src_lambda_min
REAL    :: src_lambda_max
REAL    :: targ_lambda_min
REAL    :: targ_lambda_max
REAL    :: grid_tol
REAL    :: delta_lambda

CHARACTER(LEN=*), PARAMETER :: routinename = 'lam_inclusion'
CHARACTER(LEN=80) :: cmessage

REAL, ALLOCATABLE :: lambda_box(:)
REAL, ALLOCATABLE :: phi_box(:)

LOGICAL :: l_warn_only

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)

IF (PRESENT(grid_tol_arg)) THEN
  grid_tol = grid_tol_arg
ELSE
  grid_tol = 0.0
END IF

! If we dont care about grid tolerance then only warn.
IF (grid_tol == rmdi) THEN
  l_warn_only = .TRUE.
  grid_tol = 0.0
ELSE
  l_warn_only = .FALSE.
END IF

num_pts = 2*(2*targ_halo_y+targ_rows)+2*(2*targ_halo_x+targ_row_len)
ALLOCATE(lambda_box(num_pts))
ALLOCATE(phi_box   (num_pts))

! Fill boundary box (start in SW corner and fill anti-clockwise)
! Fill south boundary
ipt = 0
DO j = 1-targ_halo_x, targ_row_len+targ_halo_x
  ipt = ipt + 1
  lambda_box(ipt) = targ_lambda(j)
  phi_box   (ipt) = targ_phi(1-targ_halo_y)
END DO
! Fill east boundary
DO j = 1-targ_halo_y, targ_rows+targ_halo_y
  ipt = ipt + 1
  lambda_box(ipt) = targ_lambda(targ_row_len+targ_halo_x)
  phi_box   (ipt) = targ_phi(j)
END DO

! Fill north boundary
DO j = targ_row_len+targ_halo_x,1-targ_halo_x,-1
  ipt = ipt + 1
  lambda_box(ipt) = targ_lambda(j)
  phi_box   (ipt) = targ_phi(targ_rows+targ_halo_y)
END DO

! Fill west boundary
DO j = targ_rows+targ_halo_y,1-targ_halo_y,-1
  ipt = ipt + 1
  lambda_box(ipt) = targ_lambda(1-targ_halo_x)
  phi_box   (ipt) = targ_phi(j)
END DO

IF (ipt /= num_pts) THEN
  WRITE(6,'(A,I5,A,I5)') "ipt = ", ipt, ", num_pts = ", num_pts
  cmessage = "Inconsistency detected filling up boundary coordinates"
  errorstatus = 1
  CALL ereport( routinename, errorstatus, cmessage )
END IF

IF (.NOT. l_same_rotation) THEN
  CALL EqToLL (                                                 &
               phi_box,  lambda_box,                            &
               phi_box,  lambda_box,                            &
               targ_pole_lat, targ_pole_lon, num_pts )
  IF (l_src_rotated) THEN

    CALL LLToEq (                                               &
      phi_box,                                                  &
      lambda_box,                                               &
      phi_box,                                                  &
      lambda_box,                                               &
      src_pole_lat,                                             &
      src_pole_lon,                                             &
      num_pts                                                   &
    )
  END IF
END IF

! Ensure we have a boundary box which is consistent.  We only care about the
! lambda due to sphere being cyclic.
lambda_box(1) = MODULO(lambda_box(1),360.)
DO ipt = 2, num_pts
  ! Calculate spacing
  delta_lambda = lambda_box(ipt)-lambda_box(ipt-1)
  ! Make sure we have the nearest point taking into account crossing meridian.
  lambda_box(ipt) = lambda_box(ipt)-NINT(delta_lambda/360.)*360.
END DO

! Make sure we are all positive so our mostly westerly point is positive.  We
! make sure the most westerly point is also positive on source grid later.
IF (MINVAL(lambda_box) < 0.0) THEN
  lambda_box(:) = lambda_box(:) + 360.
END IF

! Store source leftmost and rightmost longitudes in temporary variables for
!   clarity
src_lambda_min = MODULO(src_lambda(1-src_halo_x),360.)
src_lambda_max = src_lambda(src_row_len+src_halo_x)
DO WHILE(src_lambda_max < src_lambda_min)
  src_lambda_max = src_lambda_max + 360.
END DO
src_phi_min    = src_phi(1-src_halo_y)
src_phi_max    = src_phi(src_rows+src_halo_y)

! Now we have lambda_box with all the lambda coordinates with reference to the
! source grid.  Lets make sure all our target lambdas are to the east of the
! source grid.

! Make sure all our target box is to the east of the source box.
DO ipt = 1, num_pts
  IF (lambda_box(ipt) + grid_tol < src_lambda_min) THEN
    lambda_box(:) = lambda_box(:) + 360.
    EXIT
  END IF
END DO
! Calculate minimum lambda.
targ_lambda_min = MINVAL(lambda_box)

! Now calculate maximum lambda.
targ_lambda_max = MAXVAL(lambda_box)


! Deallocate temp arrays
DEALLOCATE(phi_box)
DEALLOCATE(lambda_box)

! ----------------------------------------
! Check target grid is inside source grid
! ----------------------------------------
! If the leftmost point in the target longitudes  is to the left of the 
!   leftmost source longitude (source - target = positive) then the target
!   is outside the source
IF ( src_lambda_min - (targ_lambda_min + grid_tol) > 0.0 ) THEN
  WRITE(cmessage,'(A,F8.3,A,F8.3)')                &
    'Western boundary at ', targ_lambda_min,       &
    ' not within source LAM with ', src_lambda_min
  
  IF (l_warn_only) THEN
    errorstatus=-4
  ELSE
    errorstatus=4
  END IF

  CALL ereport( routinename, errorstatus, cmessage )
END IF


! If the rightmost point in the target longitudes  is to the right of the 
!   rightmost source longitude (source - target = negative) then the target
!   is outside the source
IF ( src_lambda_max - (targ_lambda_max - grid_tol) < 0.0 ) THEN
  WRITE(cmessage,'(A,F8.3,A,F8.3)')                &
    'Eastern boundary at ', targ_lambda_max,       &
    ' not within source LAM with ', src_lambda_max
  
  IF (l_warn_only) THEN
    errorstatus=-4
  ELSE
    errorstatus=4
  END IF

  CALL ereport( routinename, errorstatus, cmessage )
END IF


! If the lowermost target point is below the lowermost source point the
!   target is outside the source
IF (targ_phi_min + grid_tol < src_phi_min ) THEN
  WRITE(cmessage,'(A,F8.3,A,F8.3)')            &
    'Southern boundary at ', targ_phi_min,     &
    ' not within source LAM with ', src_phi_min 
  
  IF (l_warn_only) THEN
    errorstatus=-4
  ELSE
    errorstatus=4
  END IF

  CALL ereport( routinename, errorstatus, cmessage )
END IF


! If the uppermost target point is above the uppermost source point the
!   target is outside the source
IF (targ_phi_max - grid_tol > src_phi_max ) THEN
  WRITE(Cmessage,'(A,F8.3,A,F8.3)')                   &
    'Northern boundary at ', targ_phi_max,            &
    ' not within source LAM with ', src_phi_max
  
  IF (l_warn_only) THEN
    errorstatus=-4
  ELSE
    errorstatus=4
  END IF

  CALL ereport( routinename, errorstatus, cmessage )
END IF

IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)
END SUBROUTINE lam_inclusion
END MODULE lam_inclusion_mod
