! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE NEAR_PT-------------------------------------------------
!LL
!LL  Purpose:  To produce gather indices which map each point on the
!LL            target grid onto the nearest point on the source grid.
!LL            This allows interpolation by choosing the value of the
!LL            nearest neighbour. The code uses coefficients and gather
!LL            indices calculated by subroutine H_INT_CO.
!LL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  5.1      10/04/00      South->North ordering fix. P.Selwood.
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  System component:
!LL
!LL  System task: S123
!LL
!LL  Documentation: The interpolation formulae are described in
!LL                 unified model on-line documentation paper S1.
!LL
!LL  -------------------------------------------------------------------
!*L  Arguments:---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

      SUBROUTINE NEAR_PT                                                &
     &(INDEX_B_L,INDEX_B_R,WEIGHT_T_R,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L  &
     &,POINTS,POINTS_LAMBDA_SRCE,INDEX_NEAREST)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS_LAMBDA_SRCE                                               &
                            !IN Number of lambda points on source grid
     &,POINTS                                                           &
                            !IN Total number of points on target grid
     &,INDEX_B_L(POINTS)                                                &
                            !IN  Index of bottom lefthand corner
                            !    of source gridbox
     &,INDEX_B_R(POINTS)                                                &
                            !IN  Index of bottom righthand corner
                            !    of source gridbox
     &,INDEX_NEAREST(POINTS)!OUT Index of nearest source point to
                            ! each target point

      REAL                                                              &
     & WEIGHT_T_R(POINTS)                                               &
                          !IN  Weight applied to value at top right
                          !    hand corner of source gridbox
     &,WEIGHT_B_L(POINTS)                                               &
                          !IN  Weight applied to value at bottom left
                          !    hand corner of source gridbox
     &,WEIGHT_B_R(POINTS)                                               &
                          !IN  Weight applied to value at bottom right
                          !    hand corner of source gridbox
     &,WEIGHT_T_L(POINTS) !IN  Weight applied to value at top left
                          !    hand corner of source gridbox

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & INDEX_TEMP(POINTS,4)   ! Index of 4 sourrounding source points
                              ! ordered by distance

      REAL                                                              &
     & MAX_WEIGHT(POINTS,4)   ! Linear interpolation weights ordered by
                              ! distance

!*L External subroutines called:----------------------------------------
! None
!*----------------------------------------------------------------------
! Local variables:------------------------------------------------------
      REAL TEMP
      INTEGER I,J,K,ITEMP
! ----------------------------------------------------------------------

! 1.  Accumulate source weights and indices associated with
!     each coastal point on target grid.

      DO I=1,POINTS

      MAX_WEIGHT(I,1)=WEIGHT_B_L(I)
      MAX_WEIGHT(I,2)=WEIGHT_B_R(I)
      MAX_WEIGHT(I,3)=WEIGHT_T_L(I)
      MAX_WEIGHT(I,4)=WEIGHT_T_R(I)
      INDEX_TEMP(I,1)=INDEX_B_L(I)
      INDEX_TEMP(I,2)=INDEX_B_R(I)
      INDEX_TEMP(I,3)=INDEX_B_L(I)                                      &
     &                +POINTS_LAMBDA_SRCE
      INDEX_TEMP(I,4)=INDEX_B_R(I)                                      &
     &                +POINTS_LAMBDA_SRCE
      END DO

! 2.  Sort gather indices of the 4 surrounding source
!     gridpoints according to distance from target gridpoint;
!     arranged so that nearest point comes first in list (ie K=1).

      DO K=1,3
      DO J=K+1,4
      DO I=1,POINTS
      IF(MAX_WEIGHT(I,K) <  MAX_WEIGHT(I,J))THEN
      TEMP=MAX_WEIGHT(I,K)
      MAX_WEIGHT(I,K)=MAX_WEIGHT(I,J)
      MAX_WEIGHT(I,J)=TEMP
      ITEMP=INDEX_TEMP(I,K)
      INDEX_TEMP(I,K)=INDEX_TEMP(I,J)
      INDEX_TEMP(I,J)=ITEMP
      ENDIF
      END DO
      END DO
      END DO

! 3. Assign index of nearest source point to output array

      DO I=1,POINTS
      INDEX_NEAREST(I)=INDEX_TEMP(I,1)
      END DO

      RETURN
      END SUBROUTINE NEAR_PT
