! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE COAST_AJ-----------------------------------------------
!
!
!    Purpose:
!         (i)  To produce gather indices which map each coastal point
!              on the target grid onto its nearest point on the source
!              grid. This allows correction of those surface fields
!              which are non-homogeneous across land/sea boundaries
!              after horizontal interpolation by subroutine H_INT.
!              The algorithm uses linear interpolation weights and
!              gather indices calculated by subroutine H_INT_CO.
!
!         (ii) If a land-sea mask for the target grid is not provided,
!              one is created. When a target land/sea mask is provided,
!              further index is output containing those points on the
!              target grid for which the 4 surrounding source points
!              are not of the same land/sea type as the target point.
!              These points will generally be new islands etc resolved b
!              a hires land/sea mask. They should be set to appropriate
!              values (eg climatology) as required.
!
!    Programming standard: UMDP3 v8.2
!
!    Logical component number: S122
!
!    Project task: S1
!
!    Documentation: The interpolation formulae are described in
!                   unified model on-line documentation paper S1.
!
!    ------------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Grids
MODULE coast_aj_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE coast_aj                                               &
(index_b_l,index_b_r,weight_t_r,weight_b_r,weight_t_l,weight_b_l  &
,points_lambda_srce,points_phi_srce,points,land_sea_srce          &
,land_sea_targ,index_targ,index_srce,coastal_points,mask          &
,index_targ_sea_unres,sea_points_unres                            &
,index_targ_land_unres,land_points_unres)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) :: points_lambda_srce    
                             !IN Number of lambda points on source grid
INTEGER, INTENT(IN) :: points_phi_srce 
                             !IN Number of phi points on source grid
INTEGER, INTENT(IN) :: points 
                             !IN Total number of points on target grid
INTEGER, INTENT(IN) :: index_b_l(points)
                             !IN  Index of bottom lefthand corner of source gridbox         
INTEGER, INTENT(IN) :: index_b_r(points) 
                             !IN  Index of bottom righthand corner of source gridbox       
INTEGER, INTENT(IN) :: land_sea_srce(points_lambda_srce*points_phi_srce)
                             !IN Land/sea mask on source grid 
                              
INTEGER, INTENT(INOUT) :: land_sea_targ(points)  
                             !INOUT Land/sea mask on target grid. If MASK
                             !then precalculated land/sea mask on input.                           

INTEGER, INTENT(OUT) :: coastal_points
                             !OUT Number of coastal points on target grid
INTEGER, INTENT(OUT) :: land_points_unres
                             !OUT No of unresolved land pts when MASK=T
INTEGER, INTENT(OUT) :: sea_points_unres 
                             !OUT No of unresolved sea pts when MASK=T

INTEGER, INTENT(OUT) :: index_targ(points)  
                             !OUT Index of target coastal points
INTEGER, INTENT(OUT) :: index_srce(points)
                             !OUT Index of source points mapped onto
                             !target coastal points
INTEGER, INTENT(OUT) :: index_targ_sea_unres(points)
                             !OUT Index of sea pts on target grid which a
                             !unresolved when MASK=T
INTEGER, INTENT(OUT) :: index_targ_land_unres(points)
                             !OUT Index of land pts on target grid which
                             !unresolved when MASK=T


REAL, INTENT(IN) ::  weight_t_r(points) 
                             !IN  Weight applied to value at top right
                             !    hand corner of source gridbox
REAL, INTENT(IN) ::  weight_b_l(points)
                             !IN  Weight applied to value at bottom left
                             !    hand corner of source gridbox
REAL, INTENT(IN) ::  weight_b_r(points)
                             !IN  Weight applied to value at bottom right
                             !    hand corner of source gridbox
REAL, INTENT(IN) ::  weight_t_l(points) 
                             !IN  Weight applied to value at top left
                             !    hand corner of source gridbox

LOGICAL, INTENT(IN) :: mask
                             !IN =F, then land/sea mask estimated
                             !   =T, then land/sea mask input as LAND_SEA_TARG


! Local arrays
INTEGER   ::   land_sea_temp(points)
                        ! Array used to accumulate output land/sea
INTEGER   ::   index_temp(points,4)
                        ! Index of 4 sourrounding source points
                        ! order by distance
INTEGER   ::   land_sea_coast(points) 
                        ! Mask of coastal points on target grid

REAL      ::   max_weight(points,4)   
                        ! Linear interpolation weights ordered by

LOGICAL   ::   logic_test(points)   ! Logical string used to accumulate tests


REAL      ::   temp
INTEGER   ::   i,j,k,start,itemp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!  1. Sum of land/sea mask values on source grid surrounding each point
!    on target grid. This assumes that the mask is integer 1's and 0's.

IF (lhook) CALL dr_hook('COAST_AJ',zhook_in,zhook_handle)
DO i=1,points

  land_sea_temp(i)=land_sea_srce(index_b_l(i))                    &
           +land_sea_srce(index_b_r(i))                           &
           +land_sea_srce(index_b_l(i)+points_lambda_srce)        &
           +land_sea_srce(index_b_r(i)+points_lambda_srce)

END DO

!  2. Generate coastal gather indices and new land/sea mask

! 2.1 Gather index for coastal points on target grid

DO i=1,points
  logic_test(i)=land_sea_temp(i) >  0.AND.land_sea_temp(i) <  4
END DO

coastal_points = 0
DO i=1,points
  IF(logic_test(i))THEN
    coastal_points=coastal_points + 1
    index_targ(coastal_points) = i
  END IF
END DO


! Set target land-sea mask to land if all 4 surrounding
! source points are also land.

DO i=1,points
  IF(land_sea_temp(i) == 4)land_sea_temp(i)=1
END DO


! 2.2 Gather index for points which differ between input
! land/sea mask and first estimate and are not first
! estimate coastal points.

IF(mask)THEN
  j=coastal_points

  DO i=1,points
    land_sea_coast(i)=0
  END DO
  
  DO i=1,coastal_points
    land_sea_coast(index_targ(i))=1
  END DO

  DO i=1,points
    logic_test(i)=land_sea_temp(i) /= land_sea_targ(i)            &
                  .AND.land_sea_coast(i) == 0
  END DO

  start=coastal_points
  coastal_points = 0
  DO i=1,points
    IF(logic_test(i))THEN
      coastal_points=coastal_points + 1
      index_targ(start+coastal_points) = i
    END IF
  END DO
  coastal_points=coastal_points+j

END IF

! 2.3 Accumulate source weights and indices associated with
!     each coastal point on target grid.

DO i=1,coastal_points

  max_weight(i,1)=weight_b_l(index_targ(i))
  max_weight(i,2)=weight_b_r(index_targ(i))
  max_weight(i,3)=weight_t_l(index_targ(i))
  max_weight(i,4)=weight_t_r(index_targ(i))
  index_temp(i,1)=index_b_l(index_targ(i))
  index_temp(i,2)=index_b_r(index_targ(i))
  index_temp(i,3)=index_b_l(index_targ(i))                        &
                  +points_lambda_srce
  index_temp(i,4)=index_b_r(index_targ(i))                        &
                  +points_lambda_srce
END DO

! 2.4 Sort gather indices of the 4 surrounding source
!     gridpoints according to distance from target gridpoint;
!     arranged so that nearest point comes first in list (ie K=1).

DO k=1,3
  DO j=k+1,4
    DO i=1,coastal_points
      IF(max_weight(i,k) <  max_weight(i,j))THEN
        temp=max_weight(i,k)
        max_weight(i,k)=max_weight(i,j)
        max_weight(i,j)=temp
        itemp=index_temp(i,k)
        index_temp(i,k)=index_temp(i,j)
        index_temp(i,j)=itemp
      END IF
    END DO
  END DO
END DO


! 2.5 Initialise source gather index as nearest point on source grid
!     to target coastal point.
DO i=1,coastal_points
  index_srce(i)=index_temp(i,1)
END DO

! 2.6 Select source gather index as nearest point of same land/sea type
!     when a land/sea mask has been input on target grid
IF(mask) THEN

  DO k=4,1,-1
    DO i=1,coastal_points
IF(land_sea_targ(index_targ(i)) == land_sea_srce(index_temp(i,k)))&
                  index_srce(i)=index_temp(i,k)
    END DO
  END DO

END IF

! 2.7 Update coastal values of land-sea mask

DO i=1,coastal_points
  land_sea_temp(index_targ(i))=land_sea_srce(index_srce(i))
END DO


! 2.8 Overwrite target land/sea mask with estimate or output
! indices of target land/sea points for which none of 4 surrounding
! source points are of same type.

IF(.NOT.mask)THEN
  DO i=1,points
    land_sea_targ(i)=land_sea_temp(i)
  END DO
  land_points_unres=0
  sea_points_unres=0

ELSE

  DO i=1,points
    logic_test(i)=land_sea_targ(i) == 0.AND.                      &
                 (land_sea_temp(i) /= land_sea_targ(i))
  END DO

  sea_points_unres = 0

  DO i=1,points
    IF(logic_test(i))THEN
      sea_points_unres=sea_points_unres + 1
      index_targ_sea_unres(sea_points_unres) = i
    END IF
  END DO

  DO i=1,points
    logic_test(i)=land_sea_targ(i) == 1.AND.                      &
                 (land_sea_temp(i) /= land_sea_targ(i))
  END DO

  land_points_unres = 0
  
  DO i=1,points
    IF(logic_test(i))THEN
      land_points_unres=land_points_unres + 1
      index_targ_land_unres(land_points_unres) = i
    END IF
  END DO

END IF

IF (lhook) CALL dr_hook('COAST_AJ',zhook_out,zhook_handle)
RETURN
END SUBROUTINE coast_aj
END MODULE coast_aj_mod
