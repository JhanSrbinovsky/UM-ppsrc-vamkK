! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Performs Bi-linear horizitontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_LSM(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT          &
     &,                    RMDI                                         &
     &,                    INDEX_B_L,INDEX_B_R,DATA_IN                  &
     &,                    WEIGHT_B_L,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_T_R  &
     &,                    LSMO                                         &
     &,                    DATA_OUT)


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
! This file belongs in section: Small execs
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.5   04/09/98   New routine: based on hint_bl. Checks input values
!                    for missing data indicator (RMDI). Output is
!                    RMDI at a point if any input point is RMDI
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER  ROWS_IN              !No of P rows on source grid
      INTEGER  ROW_LENGTH_IN        !No of pts per row on source grid
      INTEGER  LEN_FIELD_OUT        !No of points on target grid
      REAL     RMDI                 !real missing data indicator

!   Array  arguments with intent(in):
      INTEGER  INDEX_B_L(LEN_FIELD_OUT)
                                     !Index of bottom lefthand corner
                                     !  of source gridbox
      INTEGER  INDEX_B_R(LEN_FIELD_OUT)
                                     !Index of bottom righthand corner
                                     !  of source gridbox
      REAL     DATA_IN(ROWS_IN*ROW_LENGTH_IN)
                                      !Data before interpolation
      REAL     WEIGHT_B_L(LEN_FIELD_OUT)
                                     !Weight applied to value at bottom
                                     !lefthand corner of source gridbox
      REAL     WEIGHT_B_R(LEN_FIELD_OUT)
                                     !Weight applied to value at bottom
                                     !righthand corner of source gridbox
      REAL     WEIGHT_T_L(LEN_FIELD_OUT)
                                     !Weight applied to value at top
                                     !lefthand corner of source gridbox
      REAL     WEIGHT_T_R(LEN_FIELD_OUT)
                                     !Weight applied to value at top
                                     !righthand corner of source gridbox
      INTEGER  LSMO(LEN_FIELD_OUT)
                                     !Ocean land sea mask

!   Array  arguments with intent(out):
      REAL     DATA_OUT(LEN_FIELD_OUT) !Data after interpolation

! Local scalars:
      INTEGER      I

! Function & Subroutine calls:
!     External None

!- End of header

!     1. Carry out horizontal interpolation using equation (2.1)

      DO I=1,LEN_FIELD_OUT

      IF (    DATA_IN(INDEX_B_L(I))  ==  RMDI .OR.                      &
     &        DATA_IN(INDEX_B_R(I))  ==  RMDI .OR.                      &
     &        DATA_IN(INDEX_B_L(I)-ROW_LENGTH_IN)  ==  RMDI .OR.        &
     &        DATA_IN(INDEX_B_R(I)-ROW_LENGTH_IN)  ==  RMDI  .OR.       &
     &        LSMO(I)  ==  1 ) THEN

        DATA_OUT(I)=RMDI

      ELSE

        DATA_OUT(I)=WEIGHT_B_L(I)*DATA_IN(INDEX_B_L(I))                 &
     &             +WEIGHT_B_R(I)*DATA_IN(INDEX_B_R(I))                 &
     &             +WEIGHT_T_L(I)*DATA_IN(INDEX_B_L(I)-ROW_LENGTH_IN)   &
     &             +WEIGHT_T_R(I)*DATA_IN(INDEX_B_R(I)-ROW_LENGTH_IN)

      END IF

      END DO

      RETURN
      END SUBROUTINE H_INT_LSM

