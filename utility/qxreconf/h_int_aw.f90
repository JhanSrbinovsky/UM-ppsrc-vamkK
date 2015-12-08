! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Performs Area weighted horizitontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_AW(ROWS_IN,ROWS_OUT                              &
     &,                   ROW_LENGTH_IN,ROW_LENGTH_OUT,GLOBAL           &
     &,                   AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP           &
     &,                   AW_COLAT_T,AW_LONG_L,DATA_IN,DATA_OUT)

!LL  System component: S121
!LL
!LL  System task: S1
!LL
!LL  Purpose:
!LL
!LL  Documentation:
!LL            The interpolation formulae are described in
!LL            unified model on-line documentation paper S1.
!LL
      IMPLICIT NONE
!
! Description:
!   Carries out bi-linear horizontal interpolation using coefficients
!   and gather indices calculated in subroutine H_INT_CO
!
! Method:
!   See UMDP S1 for full desciption
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: S121
! System Task:              S1
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER      ROWS_IN          !No of rows on source grid
      INTEGER      ROWS_OUT         !No of rows on target grid
      INTEGER      ROW_LENGTH_IN    !No of pts per row on source grid
      INTEGER      ROW_LENGTH_OUT   !No of pts per row on target grid
      LOGICAL      GLOBAL           !True if global area required

!   Array  arguments with intent(in):
      INTEGER      AW_INDEX_TARG_LHS(ROW_LENGTH_OUT+1)
                                    !Index of source box overlapping
                                    !lhs of target grid-box
      INTEGER      AW_INDEX_TARG_TOP(ROWS_OUT+1)
                                    !Index of source box overlapping
                                    !top of target grid-box
      REAL         AW_COLAT_T(ROWS_OUT+1)
                                    !Colatitude of top of target grd-box
                                    ! (in units of DELTA_LAT_SRCE)
      REAL         AW_LONG_L(ROW_LENGTH_OUT+1)
                                    !Left longitude of target grid-box
                                    ! (in units of DELTA_LONG_SRCE)
      REAL         DATA_IN(ROW_LENGTH_IN*ROWS_IN)
                                    !Data before interpolation

!   Array  arguments with intent(out):
      REAL         DATA_OUT(ROW_LENGTH_OUT*ROWS_OUT)
                                    !Data after interpolation

! Local scalars:
      INTEGER      I

! Local arrays:
      REAL         BOXSUM(ROW_LENGTH_OUT,ROWS_OUT)
                                    !Sum of data on target grid

!- End of header

!     1. Calculate sum of contribution from gridboxes

! DEPENDS ON: box_sum
      CALL BOX_SUM(ROW_LENGTH_IN,ROWS_IN,ROW_LENGTH_OUT,ROWS_OUT        &
     &,            AW_LONG_L,AW_COLAT_T                                 &
     &,            AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP                  &
     &,            GLOBAL,DATA_OUT,DATA_IN)


      RETURN
      END SUBROUTINE H_INT_AW

