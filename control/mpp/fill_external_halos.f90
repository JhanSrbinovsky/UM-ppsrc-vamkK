! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine for filling in External Halos

SUBROUTINE fill_external_halos(                                   &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

USE dynamics_input_mod, ONLY: l_endgame

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

! Purpose:
!  Fills external halos (those around the edge of model domain) with
!  sensible (copy of interior points) numbers.
!  This is useful for LAM models when SWAPBOUNDS doesn't fill these
!  halos as they are normally filled with LBC data. However, not all
!  fields have LBC data applied to them, and such fields need to have
!  sensible numbers put in these external halo regions.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Arguments:

INTEGER, INTENT(IN) :: row_length ! IN: number of points on a row
                                  !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! IN: number of rows in a theta field
                                  !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! IN: number of model levels
INTEGER, INTENT(IN) :: halo_x     ! IN: size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! IN: size of halo in "j" direction

REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
                             1-halo_y:rows+halo_y,                &
                                                 levels)    
                                  ! IN/OUT : Field to have its halos updated

! Local variables

INTEGER  ::  i,j,k            ! loop indicies
INTEGER  ::  j_start,j_stop   ! loop indicies

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook('FILL_EXTERNAL_HALOS',zhook_in,zhook_handle)

! Loop limits changed to stop use of undefined data 
IF (neighbour(pwest) == nodomain .or. neighbour(peast) == nodomain) THEN
  IF (neighbour(psouth)  ==  nodomain) THEN
    j_start=1
  ELSE
    j_start=1-halo_y
  ENDIF
  IF (neighbour(pnorth)  ==  nodomain) THEN
    j_stop=rows
  ELSE
    j_stop=rows+halo_y
  ENDIF
ENDIF

IF ( l_endgame ) THEN
! parameters used pwest, peast, psouth, pnorth, nodomain
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                          &
!$OMP& SHARED(neighbour,levels,rows,row_length,halo_y,halo_x,         &
!$OMP& field,j_start,j_stop)

! Western halo region
IF (neighbour(pwest)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = j_start, j_stop
      DO i = 1-halo_x, 0
        field(i,j,k)=field(1,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PWest)  ==  NoDomain)

! Eastern halo region
IF (neighbour(peast)  ==  nodomain) THEN
 
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = j_start, j_stop
      DO i = row_length+1, row_length+halo_x
        field(i,j,k)=field(row_length,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PEast)  ==  NoDomain)

! Northern halo region
IF (neighbour(pnorth)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = rows+1, rows+halo_y
      DO i = 1-halo_x, row_length+halo_x
        field(i,j,k)=field(i,rows,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PNorth)  ==  NoDomain)

! Southern halo region
IF (neighbour(psouth)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = 1-halo_y, 0
      DO i = 1-halo_x, row_length+halo_x
        field(i,j,k)=field(i,1,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PSouth)  ==  NoDomain)

! End of OpenMP parallel region
!$OMP END PARALLEL

ELSE

! parameters used pwest, peast, psouth, pnorth, nodomain
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                          &
!$OMP& SHARED(neighbour,levels,rows,row_length,halo_y,halo_x,         &
!$OMP& field,j_start,j_stop)

! Western halo region
IF (neighbour(pwest)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = j_start, j_stop
      DO i = 1-halo_x, 0
        field(i,j,k)=field(i+halo_x,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PWest)  ==  NoDomain)

! Eastern halo region
IF (neighbour(peast)  ==  nodomain) THEN
 
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = j_start, j_stop
      DO i = row_length+1, row_length+halo_x
        field(i,j,k)=field(i-halo_x,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PEast)  ==  NoDomain)

! Northern halo region
IF (neighbour(pnorth)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = rows+1, rows+halo_y
      DO i = 1-halo_x, row_length+halo_x
        field(i,j,k)=field(i,j-halo_y,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PNorth)  ==  NoDomain)

! Southern halo region
IF (neighbour(psouth)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = 1-halo_y, 0
      DO i = 1-halo_x, row_length+halo_x
        field(i,j,k)=field(i,j+halo_y,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PSouth)  ==  NoDomain)

! End of OpenMP parallel region
!$OMP END PARALLEL

END IF  ! l_endgame

IF (lhook) CALL dr_hook('FILL_EXTERNAL_HALOS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE fill_external_halos
