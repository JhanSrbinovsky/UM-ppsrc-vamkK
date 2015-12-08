! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Thins out a field if points are below a threshold and they have a specified
! number of neighbours which are also below that threshold.
!
! Subroutine Interface:
SUBROUTINE thinfield (fieldin,fieldout,factor,ErrorStatus)

  USE IO_Mod, ONLY:         &
    PP_Header_type,         &
    PP_Field_type
  USE Err_Mod, ONLY:        &
    StatusOK
  USE FldCodes_mod, ONLY:   &
    UMSectnNo

  IMPLICIT NONE
  !
  ! Description:
  !   Routine to take an array and thin out points which are adjacent to
  !   a specified number or lower also with a non-zero value where the
  !   value is less than a specified threshold.
  !   Note that this does not deal with cyclic north-south boundaries.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: FieldCalc
  !
  ! Code description:
  !   Language: Fortran 90.
  !   This code is written to UM programming standards version 8.3.

  ! Subroutine arguments

  REAL, INTENT(IN) :: Factor(3)     ! Parameters set in namelist

  TYPE(PP_Field_type), INTENT(IN)  :: fieldin  ! Field to thin
  TYPE(PP_Field_type), INTENT(OUT) :: fieldout ! Thinned field

  INTEGER, INTENT(INOUT) :: ErrorStatus

  ! Local variables

  CHARACTER(LEN=*), PARAMETER :: RoutineName = "thinfield"
  INTEGER :: maxadjacent      ! Max Number of adjacent points to have non-zero
                              ! values for it to be thinned
  REAL    :: threshold        ! value of field below which to thin
  INTEGER :: maxqualify

  INTEGER :: i,j              ! loop counters
  INTEGER :: imask            ! denotes whether to thin a point or not
  INTEGER, ALLOCATABLE :: adjqualify(:,:)
  LOGICAL, ALLOCATABLE :: within_threshold(:,:)
  INTEGER :: numchecked       ! counters
  INTEGER :: imax, jmax       ! array bounds
  LOGICAL :: cyclic_ew        ! domain is cyclic in east-west direction

  ! End of header

! DEPENDS ON: timer
  CALL Timer( RoutineName, 3 )

  IF ( ErrorStatus /= StatusOK ) THEN
    ! Previous error - do not proceed
    GO TO 9999
  END IF


! Set parameters for thinning
  maxadjacent = INT(factor(1))    ! Maximum number of adjacent points which
                                  ! can meet the thinning criteria before the
                                  ! field isn't 'small-scale'
  threshold   = factor(2)         ! Maximum value above which the code won't
                                  ! thin
  maxqualify  = INT(factor(3))    ! If a neighbour has a adjqualify above this
                                  ! number, it won't be thinned

! Get array bounds
  imax = fieldin % Hdr % NumCols
  jmax = fieldin % Hdr % NumRows


! Set the output array
  fieldout % Hdr = fieldin % Hdr
  NULLIFY ( fieldout % RData )
  ALLOCATE( fieldout % RData(fieldout % Hdr % NumCols, &
                          fieldout % Hdr % NumRows) )

! Allocate local arrays
  ALLOCATE(adjqualify(imax,jmax),within_threshold(imax,jmax))

! Check for EW-cyclic boundaries (NS not handled)
  IF ( fieldin % Hdr % LBHem == 0 .OR. fieldin % Hdr % LBHem == 4) THEN
    cyclic_ew = .TRUE.
  ELSE
    cyclic_ew = .FALSE.
  END IF

! Generate an array of points which meet the field value threshold criterion
! Don't use WHERE as the IBM compiler doesn't optimise it very well
  DO j = 1, jmax
    DO i = 1, imax
      IF (fieldin % RData(i,j) > 0.0 .AND. fieldin % RData(i,j) < threshold)  &
      THEN
          within_threshold(i,j) = .TRUE.
      ELSE
          within_threshold(i,j) = .FALSE.
      END IF
    END DO
  END DO


! For each point, count the number of neighbouring points which meet the
! threshold criterion
  DO j = 1, jmax
    DO i = 1, imax

      adjqualify(i,j) = 0
    ! Search adjacent points

! The MODULO part deals with the wrap around the e-w direction, and is
! effectively saying i = MODULO(i-1,imax)+1 for each adjacent point
! In the first argument for MODULO to compensate for FORTRAN array indices
! beginning at 1:
!   "i" is (i+1)-1  for the point to the right
!   "i-2" is (i-1)-1  for the point to the left
! (MODULO has the correct behaviour for negative arguments)

    ! Left
      IF ( i > 1 .OR. cyclic_ew ) THEN
        IF (within_threshold(MODULO(i-2,imax)+1,j) ) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Right
      IF ( i < imax .OR. cyclic_ew ) THEN
        IF (within_threshold(MODULO(i,imax)+1,j) ) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Bottom
      IF ( j > 1 ) THEN
        IF (within_threshold(i,j-1)) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Top
      IF ( j < jmax ) THEN
        IF (within_threshold(i,j+1)) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Bottom left
      IF ( (i > 1 .OR. cyclic_ew) .AND. j > 1 ) THEN
        IF (within_threshold(MODULO(i-2,imax)+1,j-1)) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Bottom right
      IF ( (i < imax .OR. cyclic_ew) .AND. j > 1 ) THEN
        IF (within_threshold(MODULO(i,imax)+1,j-1)) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Top right
      IF ( (i < imax .OR. cyclic_ew) .AND. j < jmax ) THEN
        IF (within_threshold(MODULO(i,imax)+1,j+1) ) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF

    ! Top left
      IF ( (i > 1 .OR. cyclic_ew) .AND. j < jmax ) THEN
        IF (within_threshold(MODULO(i-2,imax)+1,j+1) ) THEN
          adjqualify(i,j) = adjqualify(i,j) + 1
        END IF
      END IF
    END DO
  END DO


! For a point to be thinned it must:
!   a) be within the threshold field value
!   b) have less than a specified number of neighbours whih meet the a)
!   c) not have any neighbours which themselves have a number of neighbours
!      meeting a) greater than another specified value
  DO j = 1, jmax
    DO i = 1, imax

      imask = 0
      IF (.NOT. within_threshold(i,j)) imask = 1
      IF (adjqualify(i,j) > maxadjacent) imask = 1

    ! Left
      IF ( i > 1 .OR. cyclic_ew ) THEN
        IF (adjqualify(MODULO(i-2,imax)+1,j) > maxqualify ) imask = 1
      END IF

    ! Right
      IF ( i < imax .OR. cyclic_ew ) THEN
        IF (adjqualify(MODULO(i,imax)+1,j)> maxqualify ) imask = 1
      END IF

    ! Bottom
      IF ( j > 1 ) THEN
        IF (adjqualify(i,j-1)> maxqualify) imask = 1
      END IF

    ! Top
      IF ( j < jmax ) THEN
        IF (adjqualify(i,j+1)> maxqualify) imask = 1
      END IF

    ! Bottom left
      IF ( (i > 1 .OR. cyclic_ew) .AND. j > 1 ) THEN
        IF (adjqualify(MODULO(i-2,imax)+1,j-1)> maxqualify) imask = 1
      END IF

    ! Bottom right
      IF ( (i < imax .OR. cyclic_ew) .AND. j > 1 ) THEN
        IF (adjqualify(MODULO(i,imax)+1,j-1)> maxqualify) imask = 1
      END IF

    ! Top right
      IF ( (i < imax .OR. cyclic_ew) .AND. j < jmax ) THEN
        IF (adjqualify(MODULO(i,imax)+1,j+1) > maxqualify) imask = 1
      END IF

    ! Top left
      IF ( (i > 1 .OR. cyclic_ew) .AND. j < jmax ) THEN
        IF (adjqualify(MODULO(i-2,imax)+1,j+1) > maxqualify) imask = 1
      END IF

! This point should be thinned if imask is still zero
      fieldout % RData(i,j) = imask * fieldin % RData(i,j)


    END DO
  END DO


  DEALLOCATE(within_threshold,adjqualify)

9999 CONTINUE

! DEPENDS ON: timer
  CALL Timer( RoutineName, 4 )

END SUBROUTINE thinfield
