! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT*************************************

MODULE s_interp_mod

!-----------------------------------------------------------------------------
! Description:
!   SCM Interpolation
!-----------------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

    USE ereport_mod, ONLY : ereport
IMPLICIT NONE

CONTAINS 

  SUBROUTINE interp1d(old_x, old_y, new_x, new_y)

    !-------------------------------------------------------------------------
    ! Method:
    !   Performs linear interpolation between single rank data arrays. Source
    !   data are sorted to increase monontonically with the independent
    !   variable before interpolation. Interpolation data points which are
    !   outside the range of the original dataset are returned with missing
    !   data indicators
    !-------------------------------------------------------------------------

    USE sort_mod
    USE scm_utils, ONLY: rmdi, imdi
    USE um_types
    USE parkind1, ONLY: jpim, jprb
    USE yomhook,  ONLY: lhook, dr_hook

    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE

    REAL, INTENT(In)    :: old_x(:)
    REAL, INTENT(In)    :: old_y(:)
    REAL, INTENT(In)    :: new_x(:)
    REAL, INTENT(InOut) :: new_y(:)



    ! Local
    CHARACTER(len=8), PARAMETER :: routinename = 'interp1d'
    CHARACTER(len=80) :: cmessage  ! Error message if icode >0
    INTEGER           :: icode    

    INTEGER, ALLOCATABLE :: lower_index(:)
    REAL,    ALLOCATABLE :: weighting(:), delta_old_y(:)
    REAL,    ALLOCATABLE :: sorted_old_x(:), sorted_old_y(:)

    INTEGER :: n_old
    INTEGER :: n_new
    INTEGER :: i,j,k

    ! Dr Hook
    !==============================
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('INTERP1D',zhook_in,zhook_handle)
    !==============================

    n_old = SIZE(old_x)
    n_new = SIZE(new_x)

    ALLOCATE( lower_index  (n_new)   )
    ALLOCATE( weighting    (n_new)   )
    ALLOCATE( delta_old_y  (n_old-1) )
    ALLOCATE( sorted_old_x (n_old)   )
    ALLOCATE( sorted_old_y (n_old)   )

    icode          = 0
    new_y(:)       = rmdi
    lower_index(:) = imdi


    IF (SIZE(old_y) /= n_old)  THEN
      icode = 510
      cmessage = '(n_old_x =/ n_old_y) for interpolation to new levels'

      
      CALL ereport (routinename, icode, cmessage)
    END IF

    IF (SIZE(new_y) /= n_new) THEN
      icode = 511
      cmessage = '(n_new_x =/ n_new_y) for interpolation to new levels'

      
      CALL ereport (routinename, icode, cmessage)
    END IF

    ! Ensure avaiable data rises monotonically in independent variable i.e. x
    sorted_old_x(:) = old_x(:)
    sorted_old_y(:) = old_y(:)

    CALL sort(sorted_old_x, sorted_old_y)

    ! Assuming x co-ord of old profile rises monotonically
    ! Calculate weightings
    DO k=1, n_old-1
      delta_old_y(k) = (sorted_old_y(k+1) - sorted_old_y(k))
    END DO

    DO k=1, n_new

      IF ( (new_x(k) <= sorted_old_x(1)) ) THEN
        ! Extrapolate below the lowest namelist level 
        ! (has to be extrapolation for exner to be sensible)
        ! Note that forcings are also extrapolated here but then 
        ! reset to be equal to the lowest namelist level value
        ! in order to conserve the vertically integrated forcing.
        ! Currently all tracer gases and aerosols are also reset to the 
        ! lowest namelist level.
        ! Initial and background profiles are NOT currently reset 
        ! and so remain extrapolated.

        weighting(k) = ( new_x(k)    - sorted_old_x(1) )                    &
                     / ( sorted_old_x(2) - sorted_old_x(1) )
        lower_index(k) = 1
!      ELSE IF ( (new_x(k) > sorted_old_x(n_old)) ) THEN
!        ! extrapolate above the highest namelist level 
!        ! (this might be more dangerous!)
!
!        weighting(k) = ( new_x(k)    - sorted_old_x(n_old-1) )              &
!                     / ( sorted_old_x(n_old) - sorted_old_x(n_old-1) )
!        lower_index(k) = n_old-1
      END IF

      DO j=1, n_old-1
        IF ((new_x(k) <= sorted_old_x(j+1)) .AND.                           &
            (new_x(k) >  sorted_old_x(j))) THEN

          weighting(k) = ( new_x(k)    - sorted_old_x(j) )                  &
                       / ( sorted_old_x(j+1) - sorted_old_x(j) )

          lower_index(k) = j

        END IF
      END DO
    END DO

    DO k=1, n_new
      IF (lower_index(k) /= imdi) THEN
        new_y(k) =  old_y(lower_index(k))                                   &
                 +  weighting(k) * delta_old_y(lower_index(k))
      END IF
    END DO

    DEALLOCATE(lower_index)
    DEALLOCATE(weighting)
    DEALLOCATE(sorted_old_x)
    DEALLOCATE(sorted_old_y)


    IF (lhook) CALL dr_hook('INTERP1D',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE interp1d

!-----------------------------------------------------------------------------
END MODULE s_interp_mod
