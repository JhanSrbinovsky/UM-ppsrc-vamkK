! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


SUBROUTINE FILL_N_DSPEC( Iter,         &  !inout
                         DSpec_Fld,    &  !inout
                         Ref_Fld,      &  !inout
                         ErrorStatus)     !inout
!
!
! Description:
!   Fills in small gaps in large CBs and smooth edges and despeckles the field.
!
! Method:
!   Fields speckled with Cbs are modified to remove groups of adjacent grid
!   points with few Cbs. Areas are eroded Iter times. This removes isolated
!   Cbs. The remaining areas are then grown to preserve the shape and size
!   of the Cb group.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!
USE IO_mod, ONLY:        &
    PP_Header_type,      &
    PP_Field_type
USE FldCodes_Mod, ONLY:  &
    ST_I_ECBB
USE Err_Mod, ONLY:        &
    StatusOK

IMPLICIT None

INTEGER, INTENT(INOUT) :: Iter      ! No. iterations/erosions

! Field to be despeckled; adjusted several
!  times during despeckling process. On
!  output contains despeckled field.
TYPE(PP_Field_type), INTENT(INOUT) :: DSpec_Fld

! Duplicate of field to be despeckled.
! Used as a reference and not changed
! during despeckling process.
! If the input field is not being compared
! to the embedded Cb mask, the final
! despeckled field is copied to this
! field at the end of the subroutine.
TYPE(PP_Field_type), INTENT(INOUT) :: Ref_Fld

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "FILL_N_DSPEC"
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

INTEGER :: i, j ,k                            ! Loop counters
INTEGER :: ii, jj                             ! Optimisation loop counters


!  Variables p contain values in a 3x3 box around the point of interest (p5)
!  and are reset to 0 if they have missing values before being used to sum
!  various parts of the 3x3 grid.
!
!  p1,p2,p3,
!  p4,p5,p6,
!  p7,p8,p9,
!
! The other values are used to contain various sums of these values.

REAL    :: p1, p2, p3, p4, p5, p6, p7, p8, p9
REAL    :: p10                   ! sum of p1 to p9
REAL    :: q1                    ! sum of top row p1+p2+p3
REAL    :: q2                    ! sum of bottom row p7+p8+p9
REAL    :: q3                    ! sum of left column p1+p4+p7
REAL    :: q4                    ! sum of right column p3+p6+p9
REAL    :: q5                    ! sum of top left hand square p1+p2+p4+p5
REAL    :: q6                    ! sum of top right hand square p2+p3+p5+p6
REAL    :: q7                    ! sum of bottom right hand square p5+p6+p8+p9
REAL    :: q8                    ! sum of bottom left hand square p4+p5+p7+p8
REAL    :: q9                    ! sum of top left hand corner p1+p2+p4
REAL    :: q10                   ! sum of top right hand corner p2+p3+p6
REAL    :: q11                   ! sum of bottom right hand corner p6+p8+p9
REAL    :: q12                   ! sum of bottom left hand corner p4+p7+p8


! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF



!==== De speckle CBs =====


!--- Stage 1: remove isolated Cbs by checking neighbouring points for Cbs ----
!
! This stage is done for all input field types.
!
! For each point in the grid, count the number of gridpoints containing Cbs
! in the 3x3 grid square. If point of interest is marked as containing Cb
! and less than 3 surrounding gridpoints contain a Cb, mark the point of
! interest as not containing a Cb.
! Repeat Iter times.
!
! Optimisation comments:
!kk Optimisation exchanged I and J loop to avoid bank conflicts
!kk In the I loop, there are dependencies in Pfields(7)%RData which
!kk prevent vectorization. Running this loop with stride 3 allows
!kk vectorization, but slightly changes the algorithm .

DO k = 1, Iter
  DO jj = 0, 2
    DO j = jj+2, DSpec_Fld % Hdr % NumRows  -1, 3
      DO ii = 0, 2
        DO i = ii+2, DSpec_Fld % Hdr % NumCols -1, 3
          p1 = DSpec_Fld % RData(i-1,j+1)
          p2 = DSpec_Fld % RData(i,j+1)
          p3 = DSpec_Fld % RData(i+1,j+1)
          p4 = DSpec_Fld % RData(i-1,j)
          p5 = DSpec_Fld % RData(i,j)
          p6 = DSpec_Fld % RData(i+1,j)
          p7 = DSpec_Fld % RData(i-1,j-1)
          p8 = DSpec_Fld % RData(i,j-1)
          p9 = DSpec_Fld % RData(i+1,j-1)
          IF (p1 == RMDI ) THEN
            p1 = 0.
          END IF
          IF (p2 == RMDI ) THEN
            p2 = 0.
          END IF
          IF (p3 == RMDI ) THEN
            p3 = 0.
          END IF
          IF (p4 == RMDI ) THEN
            p4 = 0.
          END IF
          IF (p5 == RMDI ) THEN
            p5 = 0.
          END IF
          IF (p6 == RMDI ) THEN
            p6 = 0.
          END IF
          IF (p7 == RMDI ) THEN
            p7 = 0.
          END IF
          IF (p8 == RMDI ) THEN
            p8 = 0.
          END IF
          IF (p9 == RMDI ) THEN
            p9 = 0.
          END IF
          p10 = p1+p2+p3+p4+p5+p6+p7+p8+p9
          IF  ( (p10 <= 3.) .and. (p5 == 1.) ) THEN
            DSpec_Fld % RData(i,j) = RMDI
          END IF
        END DO
      END DO
    END DO
  END DO
END DO



!---- Stage 2 - Grow areas of Cbs ----
!
! This stage is not done for fields being checked with embedded Cb mask.
!
! Loops over the field to increase the area of Cbs, and for each gridpoint,
! check 3x3 grid square surrounding it and adjust as follows:
! 1) fill middle row with Cbs if top or bottom row contains Cbs at every point
! 2) fill middle column with Cbs if left or right column contains Cbs at every
!    point
! 3) fill bottom left - top right diagonal with Cbs if top left or bottom right
!    corner contains Cbs at every point
! 4) fill top left - bottom right diagonal with Cbs if top right or bottom left
!    corner contains Cbs at every point
! 5) fill all points in 3x3 square with Cbs if top left hand, top right hand,
!    bottom left hand or bottom right hand square corner contains Cbs at every
!    point
!
! Finally, copy adjusted field DSpec_Fld into Ref_Fld.
!

IF (DSpec_Fld % Hdr % STCode /= ST_I_ECBB ) THEN
  DO k = 1,2
    DO j = 2, Ref_Fld % Hdr % NumRows  -1
      DO i = 2, Ref_Fld % Hdr % NumCols -1
        p1 = Ref_Fld % RData(i-1,j+1)
        p2 = Ref_Fld % RData(i  ,j+1)
        p3 = Ref_Fld % RData(i+1,j+1)
        p4 = Ref_Fld % RData(i-1,j  )
        p5 = Ref_Fld % RData(i  ,j  )
        p6 = Ref_Fld % RData(i+1,j  )
        p7 = Ref_Fld % RData(i-1,j-1)
        p8 = Ref_Fld % RData(i  ,j-1)
        p9 = Ref_Fld % RData(i+1,j-1)
        IF (p1 == RMDI ) THEN
          p1 = 0.
        END IF
        IF (p2 == RMDI ) THEN
          p2 = 0.
        END IF
        IF (p3 == RMDI ) THEN
          p3 = 0.
        END IF
        IF (p4 == RMDI ) THEN
          p4 = 0.
        END IF
        IF (p5 == RMDI ) THEN
          p5 = 0.
        END IF
        IF (p6 == RMDI ) THEN
          p6 = 0.
        END IF
        IF (p7 == RMDI ) THEN
          p7 = 0.
        END IF
        IF (p8 == RMDI ) THEN
          p8 = 0.
        END IF
        IF (p9 == RMDI ) THEN
          p9 = 0.
        END IF
        q1  = p1+p2+p3
        q2  = p7+p8+p9
        q3  = p1+p4+p7
        q4  = p3+p6+p9
        q5  = p1+p2+p4+p5
        q6  = p2+p3+p5+p6
        q7  = p5+p6+p8+p9
        q8  = p4+p5+p7+p8
        q9  = p1+p2+p4
        q10 = p2+p3+p6
        q11 = p6+p8+p9
        q12 = p4+p7+p8

        !Check top and bottom row
        IF  ( (q1 == 3.) .OR. (q2 == 3.) ) THEN
          DSpec_Fld % RData(i-1,j) = 1.
          DSpec_Fld % RData(i  ,j) = 1.
          DSpec_Fld % RData(i+1,j) = 1.
        END IF

        !Check left and right column
        IF  ( (q3 == 3.) .OR. (q4 == 3.) ) THEN
          DSpec_Fld % RData(i,j-1) = 1.
          DSpec_Fld % RData(i,j  ) = 1.
          DSpec_Fld % RData(i,j+1) = 1.
        END IF

        !Check top left and bottom right corner
        IF  ( (q9 == 3.) .OR. (q11 == 3.) ) THEN
          DSpec_Fld % RData(i-1,j-1) = 1.
          DSpec_Fld % RData(i  ,j  ) = 1.
          DSpec_Fld % RData(i+1,j+1) = 1.
        END IF

        !Check bottom left and top right corner
        IF  ( (q10 == 3.) .OR. (q12 == 3.) ) THEN
          DSpec_Fld % RData(i-1,j+1) = 1.
          DSpec_Fld % RData(i  ,j  ) = 1.
          DSpec_Fld % RData(i+1,j-1) = 1.
        END IF

        !Check top left square
        IF  ( (q5 == 4.)  ) THEN
          DSpec_Fld % RData(i+1,j+1) = 1.
          DSpec_Fld % RData(i+1,j  ) = 1.
          DSpec_Fld % RData(i+1,j-1) = 1.
          DSpec_Fld % RData(i  ,j-1) = 1.
          DSpec_Fld % RData(i-1,j-1) = 1.
        END IF

        !Check top right square
        IF  ( (q6 == 4.)  ) THEN
          DSpec_Fld % RData(i-1,j+1) = 1.
          DSpec_Fld % RData(i-1,j  ) = 1.
          DSpec_Fld % RData(i+1,j-1) = 1.
          DSpec_Fld % RData(i  ,j-1) = 1.
          DSpec_Fld % RData(i-1,j-1) = 1.
        END IF

        !Check bottom right square
        IF  ( (q7 == 4.)  ) THEN
          DSpec_Fld % RData(i-1,j-1) = 1.
          DSpec_Fld % RData(i-1,j  ) = 1.
          DSpec_Fld % RData(i+1,j+1) = 1.
          DSpec_Fld % RData(i  ,j+1) = 1.
          DSpec_Fld % RData(i-1,j+1) = 1.
        END IF

        !Check bottom left square
        IF  ( (q8 == 4.)  ) THEN
          DSpec_Fld % RData(i+1,j-1) = 1.
          DSpec_Fld % RData(i+1,j  ) = 1.
          DSpec_Fld % RData(i+1,j+1) = 1.
          DSpec_Fld % RData(i  ,j+1) = 1.
          DSpec_Fld % RData(i-1,j+1) = 1.
        END IF
      END DO
    END DO
    Ref_Fld % RData(:,:) = DSpec_Fld % RData(:,:)
  END DO
END IF

9999 CONTINUE

END SUBROUTINE FILL_N_DSPEC
