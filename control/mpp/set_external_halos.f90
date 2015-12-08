! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine for setting the External Halos to a provided
! constant value.

SUBROUTINE set_external_halos(                                    &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y, value)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

! Purpose:
!  Set the external halos (those around the edge of model domain) to
!  a provided constant value.
!  This is useful for LAM models because SWAPBOUNDS doesn't fill these
!  halos, they are normally filled with LBC data. However, not all
!  fields have LBC data applied to them, and such fields need to have
!  sensible numbers put in these external halo regions.
!  This subroutine is similar in principle to FILL_EXTERNAL_HALOS, but
!  a provided value rather than the adjacent point's value is used.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Arguments:

INTEGER, INTENT(IN) ::                                            &
  row_length                                                      &
                   ! number of points on a row
                   !        (not including halos)
, rows                                                            &
                   ! number of rows in a theta field
                   !        (not including halos)
, levels                                                          &
                   ! number of model levels
, halo_x                                                          &
                   ! size of halo in "i" direction
, halo_y           ! size of halo in "j" direction

REAL, INTENT(IN) ::                                               &
  value            ! value to set halos to.

REAL, INTENT(INOUT) ::                                            &
  field(1-halo_x:row_length+halo_x,                               &
        1-halo_y:rows+halo_y,                                     &
        levels)    ! Field to have its halos updated


! Local variables

INTEGER ::                                                        &
  i,j,k            ! loop indicies

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------

! Western halo region

IF (lhook) CALL dr_hook('SET_EXTERNAL_HALOS',zhook_in,zhook_handle)
IF (neighbour(pwest)  ==  nodomain) THEN
  DO k=1,levels
    DO j=1-halo_y,rows+halo_y
      DO i=1-halo_x,0
        field(i,j,k)=value
      END DO ! i
    END DO ! j
  END DO ! k
END IF ! IF (neighbour(PWest)  ==  NoDomain)

! Eastern halo region

IF (neighbour(peast)  ==  nodomain) THEN
  DO k=1,levels
    DO j=1-halo_y,rows+halo_y
      DO i=row_length+1,row_length+halo_x
        field(i,j,k)=value
      END DO ! i
    END DO ! j
  END DO ! k
END IF ! IF (neighbour(PEast)  ==  NoDomain)

! Northern halo region

IF (neighbour(pnorth)  ==  nodomain) THEN
  DO k=1,levels
    DO j=rows+1,rows+halo_y
      DO i=1-halo_x,row_length+halo_x
        field(i,j,k)=value
      END DO ! i
    END DO ! j
  END DO ! k
END IF ! IF (neighbour(PNorth)  ==  NoDomain)

! Southern halo region

IF (neighbour(psouth)  ==  nodomain) THEN
  DO k=1,levels
    DO j=1-halo_y,0
      DO i=1-halo_x,row_length+halo_x
        field(i,j,k)=value
      END DO ! i
    END DO ! j
  END DO ! k
END IF ! IF (neighbour(PSouth)  ==  NoDomain)

IF (lhook) CALL dr_hook('SET_EXTERNAL_HALOS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_external_halos

